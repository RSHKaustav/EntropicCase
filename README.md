# EntropicCase
A Python program for simulating a heat pump using hourly operational data. The program combines **TESPy** for thermodynamic modelling and **CoolProp** for fast timeseries simulation across all the timesteps

---

## CoolProp and TESPy

Running TESPy across the entire would take upwards of 20 minutes, hence TESPy is run once at a representative data point to extract heat exchanger approach temperatures. These approach temperatures were used to calibrate the hourly temperature input, that we used in CoolProp, which resulted in a faster and accurate data. 

---

## Configuration 

|Parameter | Value |
|---|---|
| Heat source fluid | Water |
| Heat source flow rate | 2 kg/s |
| Heat source outlet temperature | 20 °C (fixed) |
| Heat sink inlet temperature | 60 °C (fixed) |
| Heat sink outlet temperature | 90 °C (fixed) |
| Compressor isentropic efficiency | 0.85 |
| Default refrigerant | R134a |
| HX pressure ratio (both sides) | 0.98 |

---

## Usage
You will be prompted to enter a refrigerant:

'''
Enter refrigerant (e.g. R134a, R410A, R290, R1234yf):
'''
---

## Supported Refrigerants
| Refrigerant | Initial p_low guess |
|---|---|
| `R134a` | 4.5 bar |
| `R410A` | 9.0 bar |
| `R290` | 6.5 bar |
| `R1234yf` | 4.8 bar |
| `R717` | 6.0 bar |
| `R744` | 50.0 bar |

Other refrigerants are also supported (as long as they are available in CoolProp), but their inital p_low would be assumed as 5.0 bar.

## Thermodynamic Model
 
### Design point (TESPy)
 
TESPy builds a full component network:
 
```
Heat Source (water) ──→ EVAPORATOR ──→ Heat Source out
                              ↑
                        Refrigerant loop:
                        CycleCloser → Evaporator → Compressor → Condenser → Valve → CycleCloser
                              ↓
Heat Sink (water) ────→ CONDENSER ──→ Heat Sink out
```
 
The design solve produces the real refrigerant saturation temperatures at the evaporator and condenser. The difference between these and the secondary fluid temperatures are the **approach temperatures**:
 
```
ΔT_evap_approach = T_source_out  - T_ref_evap   (≈ 8.14 K for R134a)
ΔT_cond_approach = T_sink_out    - T_ref_cond   (≈ 4.69 K for R134a)
```
 
### Time-series (CoolProp)
 
For each hourly timestep, the refrigerant saturation temperatures are calculated using the calibrated approach temperatures:
 
```
T_evap = min(T_source_in, T_source_out) - ΔT_evap_approach
T_cond = max(T_sink_in,   T_sink_out)   - ΔT_cond_approach
```
 
The vapour-compression cycle is then solved analytically through four state points:
 
| State | Location | Condition |
|---|---|---|
| h1 | Evaporator exit | Saturated vapour (x = 1) |
| h2 | Compressor exit | Superheated vapour (actual, accounting for η_s) |
| h3 | Condenser exit | Saturated liquid (x = 0) |
| h4 | Valve exit | Isenthalpic expansion (h4 = h3) |

### Skipped Timesteps

Any timestep where Q_cond < 10 kW, were skipped as they can assumed to be effectively off. 

## Outputs
- **heat_pump_performance.csv** 
- **heat_pump_performance.png**

## Limitations 
- **Off-Design Mode** - TESPy wasn't able to converge for the given off-design conditions, as the resulting pressure was too high. 
- **Single Design Point** - TESPy is only solved at single point (57 C Heat Source Inlet)

