import warnings
from weakref import ref
warnings.filterwarnings("ignore")


import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
from tespy.networks import Network
from tespy.components import CycleCloser, Compressor, Valve, HeatExchanger, Source, Sink
from tespy.connections import Connection
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec


class HeatPump:

    PRESSURE_GUESSES = {
    "R134a":   4.5,
    "R410A":   9.0,
    "R290":    6.5,   
    "R1234yf": 4.8,
    "R717":    6.0,   
    "R744":   50.0,   
    }

    def __init__(self, refrigerant: str = "R134a", eta_s: float = 0.85):
        
        self._validate_refrigerant(refrigerant)
        self.REF = refrigerant
        self.eta_s = eta_s
        self.cp_water = 4.182

        self.timestamp    = None
        self.T_src_in     = None
        self.T_src_out    = None
        self.m_src        = None
        self.T_snk_in     = None
        self.T_snk_out    = None
        self.Q_cond_kW    = None

        self.design_point   = None
        self.offdesign_point = None
        self.dt_evap_app    = None        
        self.dt_cond_app    = None        
        self.Q_cond_design  = None         

        self.df = None
        self.valid  = None

    def load_data(self, filepath: str):
        x1 = pd.read_excel(filepath, sheet_name=None)
        self.timestamp = pd.to_datetime(x1["Heat source"]["start measurement"])
        self.T_src_in = x1["Heat source"]["T_in[degC"].values
        self.T_src_out = float(x1["Heat source"]["T_out[degC]"].iloc[0])
        self.m_src      = float(x1["Heat source"]["flow[kg/s]"].iloc[0])     
        self.T_snk_in   = float(x1["Heat sink"]["T_in[degC"].iloc[0])      
        self.T_snk_out  = float(x1["Heat sink"]["T_out[degC]"].iloc[0])    
        self.Q_cond_kW  = x1["Heat sink"]["Energy[kWh]"].values

    def caliberate(self, T_hs_in_assumed: float = 57.0, Q_cond_initial: float = 503.22):
        
        m_hk = Q_cond_initial/(self.cp_water * (self.T_snk_out - self.T_snk_in))
        self.design_point = self._tespy_solve(
            T_hs_in = T_hs_in_assumed, 
            T_hs_out = self.T_src_out, 
            m_hs = self.m_src, 
            T_hk_in = self.T_snk_in, 
            T_hk_out = self.T_snk_out, 
            m_hk = m_hk
        )

        self.dt_evap_app = self.T_src_out - self.design_point['T_ref_evap_C']
        self.dt_cond_app = self.T_snk_out - self.design_point['T_ref_cond_C']
        self.Q_cond_design = self.design_point["Q_cond_kW"]

        dp = self.design_point
        print(f"  COP          = {dp['COP']:.3f}")
        print(f"  Q_cond       = {dp['Q_cond_kW']:.1f} kW")
        print(f"  Q_evap       = {dp['Q_evap_kW']:.1f} kW")
        print(f"  P_comp       = {dp['P_comp_kW']:.1f} kW")
        print(f"  T_ref_evap   = {dp['T_ref_evap_C']:.2f} °C")
        print(f"  T_ref_cond   = {dp['T_ref_cond_C']:.2f} °C")

        # Off-design validation
        print("\nRunning TESPy off-design validation (52 °C source, 950 kW) ...")
        m_hk_od = 950.0 / (self.cp_water * (self.T_snk_out - self.T_snk_in))
        self.offdesign_point = self._tespy_solve(
            T_hs_in  = 52.0,
            T_hs_out = self.T_src_out,
            m_hs     = self.m_src,
            T_hk_in  = self.T_snk_in,
            T_hk_out = self.T_snk_out,
            m_hk     = m_hk_od,
        )
        od = self.offdesign_point
        print(f"  COP (off-design) = {od['COP']:.3f}")
        print(f"  P_comp           = {od['P_comp_kW']:.1f} kW")
        print(f"  Q_evap           = {od['Q_evap_kW']:.1f} kW")

    def run_timeseries(self, T_hs_in_assumed: float = 57.0):
        N = len(self.timestamp)
        keys = ["COP", "Q_cond_kW", "Q_evap_kW", "P_comp_kW", "T_ref_evap_C", "T_ref_cond_C", "p_low_bar", "p_high_bar", "m_ref_kgs", "COP_carnot", "pratio", "m_hk_kgs"]

        nan__row = {key: np.nan for key in keys}

        records = []
        skip = 0
        for i in range(N):
            r = self._coolprop_calculation(float(self.T_src_in[i]), self.T_src_out, self.Q_cond_kW[i], self.T_snk_in, self.T_snk_out, self.eta_s)
            if r is None:
                skip += 1
                records.append(nan__row.copy())
            else:
                records.append(r)

        print(f"Skipped {skip} rows due to insufficient Q_cond_kW values.")

        self.df = pd.DataFrame(records, index=self.timestamp)
        self.df["T_src_in_C"]       = self.T_src_in
        self.df["T_src_out_C"]      = self.T_src_out
        self.df["Q_cond_demand_kW"] = self.Q_cond_kW
        self.df["is_offdesign"]     = (
            (self.df["T_src_in_C"] < T_hs_in_assumed - 5) |
            (self.df["Q_cond_kW"]  < 0.9 * self.Q_cond_design)
        )

        self.valid = self.df.dropna(subset=["COP"])

        total_E = self.valid["P_comp_kW"].sum()
        total_Q = self.valid["Q_cond_kW"].sum()
        print(f"\n  Annual SCOP      : {total_Q/total_E:.3f}")
        print(f"  Mean COP         : {self.valid['COP'].mean():.3f}")
        print(f"  Total heat deliv.: {total_Q/1e3:.0f} MWh")
        print(f"  Total elec. input: {total_E/1e3:.0f} MWh")

    def plot(self, output_path:str = "heat_pump_performance.png"):
        
        fig = self._build_plot()
        fig.savefig(output_path, dpi=155, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f"  Figure saved → {output_path}")
    
    def save_results(self, output_path: str = "heat_pump_performance.csv"):
        self.df.to_csv(output_path)
        print(f"  Performance data saved → {output_path}")


    @staticmethod

    def _validate_refrigerant(refrigerant: str):
        if refrigerant not in CP.FluidsList():
            raise ValueError(f"'{refrigerant}' not available in CoolProp. Please choose from: {CP.fluidsList()}")

    def _tespy_solve(self, T_hs_in, T_hs_out, m_hs, T_hk_in, T_hk_out, m_hk):
        ref = self.REF
        nw = Network(fluids=[ref], T_unit='C', p_unit='bar', h_unit='kJ/kg', m_unit='kg/s')
        
        cc = CycleCloser('cc')
        comp = Compressor('comp')
        cond = HeatExchanger('cond')
        valv = Valve('valv')
        evap = HeatExchanger('evap')
        si = Source('si')  #Heat Source in
        so = Sink('so') #Heat Source out    
        ki = Source('ki') #Heat Sink in
        ko = Sink('ko') #Heat Sink out

        c1=Connection(cc,"out1",evap,"in2",label="1")
        c2=Connection(evap,"out2",comp,"in1",label="2")
        c3=Connection(comp,"out1",cond,"in1",label="3")
        c4=Connection(cond,"out1",valv,"in1",label="4")
        c5=Connection(valv,"out1",cc,"in1",label="5")
        c6=Connection(si,"out1",evap,"in1",label="6")
        c7=Connection(evap,"out1",so,"in1",label="7")
        c8=Connection(ki,"out1",cond,"in2",label="8")
        c9=Connection(cond,"out2",ko,"in1",label="9")
        nw.add_conns(c1,c2,c3,c4,c5,c6,c7,c8,c9)

        comp.set_attr(eta_s=self.eta_s)
        cond.set_attr(pr1=0.98, pr2=0.98)
        evap.set_attr(pr1=0.98, pr2=0.98)
        c8.set_attr(fluid={self.REF:0,"water":1}, T=T_hk_in, p=2, m=m_hk)
        c9.set_attr(T=T_hk_out)
        c6.set_attr(fluid={self.REF:0,"water":1}, T=T_hs_in, p=2, m=m_hs)
        c7.set_attr(T=T_hs_out)
        p_guess = self.PRESSURE_GUESSES.get(ref, 5.0)
        c1.set_attr(fluid={ref: 1, "water": 0}, p=p_guess)
        c2.set_attr(x=1)   # saturated vapour at evaporator exit
        c4.set_attr(x=0)   # saturated liquid at condenser exit

        nw.set_attr(iterinfo=False)
        nw.solve(mode="design")
        return {
            "COP"          : abs(cond.Q.val) / comp.P.val,
            "Q_cond_kW"    : abs(cond.Q.val) / 1e3,
            "Q_evap_kW"    : abs(evap.Q.val) / 1e3,
            "P_comp_kW"    : comp.P.val / 1e3,
            "T_ref_evap_C" : c2.T.val,
            "T_ref_cond_C" : c4.T.val,
            "p_low_bar"    : c1.p.val,
            "p_high_bar"   : c3.p.val,
            "m_ref_kgs"    : c1.m.val,
        }
    
    def _coolprop_calculation(self, T_hs_in, T_hs_out, Q_cond_kW,T_hk_in, T_hk_out, eta_s = 0.85):

        ref = self.REF

        if Q_cond_kW < 10 or np.isnan(Q_cond_kW):
            return None
    
        T_evap = min(T_hs_in, self.T_src_out) - self.dt_evap_app
        T_cond = max(self.T_snk_out, self.T_snk_in) + self.dt_cond_app
        T_evap_K = T_evap + 273.15
        T_cond_K = T_cond + 273.15

        p_evap = CP.PropsSI('P', 'T', T_evap_K, 'Q', 1, ref) 
        p_cond = CP.PropsSI('P', 'T', T_cond_K, 'Q', 0, ref)
        h1 = CP.PropsSI('H', 'T', T_evap_K, 'Q', 1, ref)
        s1 = CP.PropsSI('S', 'T', T_evap_K, 'Q', 1, ref)
        h2s = CP.PropsSI('H', 'P', p_cond, 'S', s1, ref)
        h2 = h1 + (h2s - h1) / eta_s
        T2= CP.PropsSI('T', 'P', p_cond, 'H', h2, ref)
        h3 = CP.PropsSI('H', 'T', T_cond_K, 'Q', 0, ref)

        m_ref = Q_cond_kW * 1e3 / (h2 - h3)
        Q_evap = m_ref * (h1 - h3) / 1e3
        W_comp = m_ref * (h2 - h1) / 1e3
        COP = Q_cond_kW / W_comp

        COP_carnot = T_cond_K / (T_cond_K - T_evap_K)
        pratio = p_cond / p_evap
        m_hk = Q_cond_kW / (self.cp_water * (T_hk_out - T_hk_in))

        return {
            "COP"          : COP,
            "Q_cond_kW"    : Q_cond_kW,
            "Q_evap_kW"    : Q_evap,
            "P_comp_kW"    : W_comp,
            "T_ref_evap_C" : T_evap,
            "T_ref_cond_C" : T_cond,
            "p_low_bar"    : p_evap / 1e5 ,  # Convert from Pa to bar
            "p_high_bar"   : p_cond / 1e5 ,   # Convert from Pa to bar
            "m_ref_kgs"    : m_ref,
            "COP_carnot"   : COP_carnot,
            "pratio"       : pratio,
            "m_hk_kgs"     : m_hk
        }
    
    def _monthly_agg(self):
        monthly = self.valid.resample("ME").agg(
            COP_mean    = ("COP",        "mean"),
            COP_min     = ("COP",        "min"),
            COP_max     = ("COP",        "max"),
            Pcomp_mean  = ("P_comp_kW",  "mean"),
            Qcond_mean  = ("Q_cond_kW",  "mean"),
            Qevap_mean  = ("Q_evap_kW",  "mean"),
        )
        mon_energy = self.valid.resample("ME").agg(
            Q_cond = ("Q_cond_kW", "sum"),
            Q_evap = ("Q_evap_kW", "sum"),
            E_comp = ("P_comp_kW", "sum"),
        ) / 1e3   # → MWh
        return monthly, mon_energy
    
    def _build_plot(self):
        valid = self.valid
        df = self.df
        dp = self.design_point
        ts = valid.index
        monthly, mon_energy = self._monthly_agg()
        mon_lbl = monthly.index.strftime("%b")

        total_E = valid["P_comp_kW"].sum()
        total_Q = valid["Q_cond_kW"].sum()
        SCOP    = total_Q / total_E
        od_hrs  = int(df["is_offdesign"].sum())
        N       = len(df)

        C = dict(
            blue="#1A5276", orange="#D35400", green="#1E8449", red="#C0392B",
            purple="#6C3483", teal="#117A8B", amber="#F39C12", gray="#808B96",
            bg="#F8F9FA", panel="#FFFFFF", text="#1C2833",
        )

        def ax_style(ax, title, ylabel=None):
            ax.set_facecolor(C["panel"])
            ax.set_title(title, fontsize=10, fontweight="bold",
                         color=C["text"], pad=5, loc="left")
            if ylabel:
                ax.set_ylabel(ylabel, fontsize=8.5, color=C["text"])
            ax.tick_params(labelsize=8, colors=C["text"])
            ax.spines[["top", "right"]].set_visible(False)
            for sp in ["left", "bottom"]:
                ax.spines[sp].set_color("#CCCCCC")
            ax.grid(axis="y", ls="--", alpha=0.35, color="#AAAAAA")
            ax.grid(axis="x", ls=":",  alpha=0.2,  color="#AAAAAA")

        def mfmt(ax):
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
            ax.tick_params(axis="x", labelsize=8)

        fig = plt.figure(figsize=(22, 25), facecolor=C["bg"])
        fig.suptitle(
            f"Heat Pump System — Full Annual Performance Analysis  (2024)\n"
            f"{self.REF}  |  Compressor η_s = {self.eta_s}  |  "
            f"Source: water {self.m_src} kg/s  |  Sink: {self.T_snk_in:.0f} → {self.T_snk_out:.0f} °C",
            fontsize=15, fontweight="bold", color=C["text"], y=0.998, va="top"
        )
        gs = GridSpec(4, 2, figure=fig, hspace=0.52, wspace=0.28,
                      left=0.07, right=0.97, top=0.965, bottom=0.03)
        
        ax1 = fig.add_subplot(gs[0, :])
        ax_style(ax1, "1  Coefficient of Performance (COP)", "COP  [–]")
        ax1.fill_between(ts, valid["COP"], alpha=0.12, color=C["blue"])
        ax1.plot(ts, valid["COP"],        color=C["blue"],  lw=0.55, label="Actual COP")
        ax1.plot(ts, valid["COP_carnot"], color=C["amber"], lw=0.5, alpha=0.65,
                 ls="--", label="Carnot COP (upper bound)")
        ax1.axhline(valid["COP"].mean(), color=C["blue"], lw=1.3, ls=":",
                    label=f"Annual mean COP = {valid['COP'].mean():.2f}")
        ax1.axhline(dp["COP"], color=C["green"], lw=1.0, ls="-.",
                    label=f"TESPy design COP = {dp['COP']:.2f}")
        ax1.fill_between(df.index, 0, 25, where=df["is_offdesign"].fillna(False),
                         alpha=0.07, color=C["red"], label="Off-design periods")
        ax1.set_ylim(0, valid["COP_carnot"].quantile(0.99) * 1.06)
        ax1.legend(fontsize=7.5, loc="upper right", ncol=3, framealpha=0.9)
        mfmt(ax1)

        # Panel 2 — Compressor power
        ax2 = fig.add_subplot(gs[1, 0])
        ax_style(ax2, "2  Compressor Power Consumption", "Power  [kW]")
        ax2.fill_between(ts, valid["P_comp_kW"], alpha=0.12, color=C["red"])
        ax2.plot(ts, valid["P_comp_kW"], color=C["red"], lw=0.55)
        ax2.axhline(valid["P_comp_kW"].mean(), color=C["red"], lw=1.2, ls=":",
                    label=f"Mean: {valid['P_comp_kW'].mean():.0f} kW")
        ax2.axhline(dp["P_comp_kW"], color=C["gray"], lw=1.0, ls="-.",
                    label=f"TESPy design: {dp['P_comp_kW']:.0f} kW")
        ax2.legend(fontsize=7.5); mfmt(ax2)

        # Panel 3 — Heat transfer rates
        ax3 = fig.add_subplot(gs[1, 1])
        ax_style(ax3, "3  Heat Transfer Rates", "Heat rate  [kW]")
        ax3.plot(ts, valid["Q_cond_kW"], color=C["orange"], lw=0.55, label="Q̇ condenser")
        ax3.plot(ts, valid["Q_evap_kW"], color=C["teal"],   lw=0.55, label="Q̇ evaporator")
        ax3.axhline(dp["Q_cond_kW"], color=C["orange"], lw=0.9, ls="-.",
                    label=f"Design Q̇_cond: {dp['Q_cond_kW']:.0f} kW")
        ax3.axhline(dp["Q_evap_kW"], color=C["teal"], lw=0.9, ls="-.",
                    label=f"Design Q̇_evap: {dp['Q_evap_kW']:.0f} kW")
        ax3.legend(fontsize=7.5, ncol=2); mfmt(ax3)

        # Panel 4 — Refrigerant temperatures
        ax4 = fig.add_subplot(gs[2, 0])
        ax_style(ax4, "4  Refrigerant Saturation Temperatures", "Temperature  [°C]")
        ax4.plot(ts, valid["T_ref_cond_C"], color=C["orange"], lw=0.55, label="T_cond (ref.)")
        ax4.plot(ts, valid["T_ref_evap_C"], color=C["teal"],   lw=0.55, label="T_evap (ref.)")
        ax4.plot(ts, valid["T_src_in_C"],   color=C["green"],  lw=0.5, alpha=0.7,
                 ls="--", label="T_source in")
        ax4.legend(fontsize=7.5, ncol=2); mfmt(ax4)

        # Panel 5 — Pressures
        ax5 = fig.add_subplot(gs[2, 1])
        ax_style(ax5, "5  Refrigerant Pressures & Pressure Ratio", "Pressure  [bar]")
        ax5.plot(ts, valid["p_high_bar"], color=C["red"],  lw=0.55, label="p_high (cond.)")
        ax5.plot(ts, valid["p_low_bar"],  color=C["teal"], lw=0.55, label="p_low (evap.)")
        ax5r = ax5.twinx()
        ax5r.plot(ts, valid["pratio"], color=C["purple"], lw=0.55, alpha=0.7,
                  label="Pressure ratio")
        ax5r.set_ylabel("Pressure ratio  [–]", fontsize=8, color=C["purple"])
        ax5r.tick_params(colors=C["purple"], labelsize=7)
        ax5r.spines["right"].set_color(C["purple"])
        h1, l1 = ax5.get_legend_handles_labels()
        h2, l2 = ax5r.get_legend_handles_labels()
        ax5.legend(h1+h2, l1+l2, fontsize=7.5); mfmt(ax5)

        # Panel 6: Monthly Average COP
        ax6 = fig.add_subplot(gs[3, 0])
        ax_style(ax6, "6  Monthly Average COP (with min/max range)", "Mean COP  [–]")
        x    = np.arange(len(monthly))
        bars = ax6.bar(x, monthly["COP_mean"], color=C["blue"], alpha=0.75,
                       edgecolor="white", width=0.6, zorder=3, label="Mean COP")
        ax6.errorbar(x, monthly["COP_mean"],
                     yerr=[monthly["COP_mean"] - monthly["COP_min"],
                           monthly["COP_max"] - monthly["COP_mean"]],
                     fmt="none", ecolor=C["gray"], elinewidth=1.2,
                     capsize=4, zorder=4, label="Min / Max range")
        ax6.set_xticks(x); ax6.set_xticklabels(mon_lbl, fontsize=8)
        for bar, v in zip(bars, monthly["COP_mean"]):
            ax6.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.04,
                     f"{v:.2f}", ha="center", va="bottom", fontsize=6.5, color=C["text"])
        ax6.axhline(SCOP, color=C["orange"], lw=1.2, ls=":",
                    label=f"Annual SCOP = {SCOP:.2f}")
        ax6.legend(fontsize=7.5)

        
        fig.text(0.5, 0.005,
                 f"Thermodynamic model: vapour-compression cycle  |  Refrigerant: {self.REF}  |  "
                 "Design point: TESPy (rigorous)  |  Time-series: CoolProp (calibrated pinch approach)",
                 ha="center", fontsize=7.5, color=C["gray"], style="italic")

        return fig


if __name__ == "__main__":
    ref = input("Enter refrigerant (e.g. R134a, R410A, R290, R1234yf): ").strip() or "R134a"

    hp = HeatPump(refrigerant=ref, eta_s=0.85)
    hp.load_data("HP_case_data.xlsx")
    hp.caliberate()
    hp.run_timeseries()
    hp.save_results("heat_pump_performance.csv")
    hp.plot("heat_pump_performance.png")





