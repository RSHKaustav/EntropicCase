import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
from tespy.networks import Network
from tespy.components import CycleCloser, Compressor, Valve, HeatExchanger, Source, Sink
from tespy.connections import Connection
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec

os.makedirs('results', exist_ok=True)

# Load data
x1 = pd.read_excel("HP_case_data.xlsx", sheet_name=None)
src_df = x1["Heat source"].copy()
snk_df = x1["Heat sink"].copy()

timestamp = pd.to_datetime(src_df["start measurement"])
T_src_in = src_df["T_in[degC"].values
T_SRC_OUT = float(src_df["T_out[degC]"].iloc[0])
M_SRC      = float(src_df["flow[kg/s]"].iloc[0])     
T_SNK_IN   = float(snk_df["T_in[degC"].iloc[0])      
T_SNK_OUT  = float(snk_df["T_out[degC]"].iloc[0])    
Q_COND_kW  = snk_df["Energy[kWh]"].values  

ETA_S = 0.85

REF = input("Enter the working fluid (e.g., R134a, R410A, etc.): ")

valid_fluids = CP.FluidsList()
if REF not in valid_fluids:
    raise ValueError(f"'{REF}' not available. Please try again")

print(f"  Using refrigerant: {REF}")

CP_Water = 4.182

# TESPy Solver
T_HS_IN_Assumed = 57.0
Q_COND_INITIAL = 503.22

def tespy_solve(T_hs_in, T_hs_out, m_hs, T_hk_in, T_hk_out, m_hk, eta_s = 0.85):
    nw = Network(fluids=[REF], T_unit='C', p_unit='bar', h_unit='kJ/kg', m_unit='kg/s')
    
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

    comp.set_attr(eta_s=eta_s)
    cond.set_attr(pr1=0.98, pr2=0.98)
    evap.set_attr(pr1=0.98, pr2=0.98)
    c8.set_attr(fluid={REF:0,"water":1}, T=T_hk_in, p=2, m=m_hk)
    c9.set_attr(T=T_hk_out)
    c6.set_attr(fluid={REF:0,"water":1}, T=T_hs_in, p=2, m=m_hs)
    c7.set_attr(T=T_hs_out)

    PRESSURE_GUESSES = {
    "R134a":   4.5,
    "R410A":   9.0,
    "R290":    6.5,   
    "R1234yf": 4.8,
    "R717":    6.0,   
    "R744":   50.0,   
    }
    p_guess = PRESSURE_GUESSES.get(REF, 5.0)  #Assumes pressure to be 5, if given refrigerant is not in the dictionary
    c1.set_attr(fluid={REF: 1, "water": 0}, p=p_guess)
    c2.set_attr(x=1)
    c4.set_attr(x=0)
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

m_hk_design = Q_COND_INITIAL/(CP_Water * (T_SNK_OUT - T_SNK_IN))
dp = tespy_solve(T_HS_IN_Assumed, T_SRC_OUT, M_SRC, T_SNK_IN, T_SNK_OUT, m_hk_design, ETA_S)

print(f"COP    = {dp['COP']:.2f}")
print(f"Q_cond = {dp['Q_cond_kW']:.2f} kW")
print(f"Q_evap = {dp['Q_evap_kW']:.2f} kW")
print(f"P_comp = {dp['P_comp_kW']:.2f} kW")
print(f"T_ref_evap = {dp['T_ref_evap_C']:.2f} °C")
print(f"T_ref_cond = {dp['T_ref_cond_C']:.2f} °C")
print(f"p_low = {dp['p_low_bar']:.2f} bar")
print(f"p_high = {dp['p_high_bar']:.2f} bar")

DeltaT_evap = T_SRC_OUT - dp['T_ref_evap_C']
DeltaT_cond = T_SNK_OUT - dp['T_ref_cond_C']

# Coolprop Calculation
def coolprop_calculation(T_hs_in, T_hs_out, Q_cond_kW,T_hk_in, T_hk_out, eta_s = 0.85):

    if Q_cond_kW < 10 or np.isnan(Q_cond_kW):
        return None
    
    T_evap = min(T_hs_out, T_hs_in) - DeltaT_evap
    T_cond = max(T_hk_out, T_hk_in) + DeltaT_cond
    T_evap_K = T_evap + 273.15
    T_cond_K = T_cond + 273.15

    p_evap = CP.PropsSI('P', 'T', T_evap_K, 'Q', 1, REF) 
    p_cond = CP.PropsSI('P', 'T', T_cond_K, 'Q', 0, REF)
    h1 = CP.PropsSI('H', 'T', T_evap_K, 'Q', 1, REF)
    s1 = CP.PropsSI('S', 'T', T_evap_K, 'Q', 1, REF)
    h2s = CP.PropsSI('H', 'P', p_cond, 'S', s1, REF)
    h2 = h1 + (h2s - h1) / ETA_S
    T2= CP.PropsSI('T', 'P', p_cond, 'H', h2, REF)
    h3 = CP.PropsSI('H', 'T', T_cond_K, 'Q', 0, REF)

    m_ref = Q_cond_kW * 1e3 / (h2 - h3)
    Q_evap = m_ref * (h1 - h3) / 1e3
    W_comp = m_ref * (h2 - h1) / 1e3
    COP = Q_cond_kW / W_comp

    COP_carnot = T_cond_K / (T_cond_K - T_evap_K)
    pratio = p_cond / p_evap
    m_hk = Q_cond_kW / (CP_Water * (T_hk_out - T_hk_in))

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

N = len(timestamp)
keys = ["COP", "Q_cond_kW", "Q_evap_kW", "P_comp_kW", "T_ref_evap_C", "T_ref_cond_C", "p_low_bar", "p_high_bar", "m_ref_kgs", "COP_carnot", "pratio", "m_hk_kgs"]

nan__row = {key: np.nan for key in keys}

records = []
skip = 0
for i in range(N):
    r = coolprop_calculation(float(T_src_in[i]), T_SRC_OUT, Q_COND_kW[i], T_SNK_IN, T_SNK_OUT, ETA_S)
    if r is None:
        skip += 1
        records.append(nan__row.copy())
    else:
        records.append(r)

print(f"Skipped {skip} rows due to insufficient Q_cond_kW values.")

#Saving results to dataframe 

df = pd.DataFrame(records, index=timestamp)
df["T_src_in_C"] = T_src_in
df["T_src_out_C"] = T_SRC_OUT
df["Q_cond_kW"] = Q_COND_kW

df.to_csv("results/heat_pump_performance.csv")

valid = df.dropna(subset = ["COP"])

monthly = valid.resample("ME").agg(
    COP_mean = ("COP", "mean"),
    COP_min = ("COP", "min"),
    COP_max = ("COP", "max"),
    Pcomp_mean = ("P_comp_kW", "mean"),
    Pcomp_max = ("P_comp_kW", "max"),
    Qcond_mean = ("Q_cond_kW", "mean"),
    Qevap_mean = ("Q_evap_kW", "mean"),

)

mon_lbl = monthly.index.strftime("%b")

mon_energy = valid.resample("ME").agg(
    Q_cond = ("Q_cond_kW", "sum"),
    Q_evap = ("Q_evap_kW", "sum"),
    E_comp = ("P_comp_kW", "sum"),
) / 1e3

total_E_comp = valid["P_comp_kW"].sum()
total_Q_cond = valid["Q_cond_kW"].sum()
total_Q_evap = valid["Q_evap_kW"].sum()
SCOP         = total_Q_cond / total_E_comp

print(f"\n  Annual SCOP        : {SCOP:.3f}")
print(f"  Mean COP           : {valid['COP'].mean():.3f}")
print(f"  Total heat deliv.  : {total_Q_cond/1e3:.0f} MWh")
print(f"  Total elec. input  : {total_E_comp/1e3:.0f} MWh")




