IRR=0.08 #annual basis
life=20 #plant lifetime years
inflation=0.038 #frac, average U.S. inflation rate from 1960 to 2021

EP = 75 # Electricity Price
Cl2P = 150 # Chlorine Price

prod_tpd=100 # tonne per day or ~40,000 tpy
employees = 44 # Employees

salary = 50000 # Assumed salary US Chemical Plant Operator
labor = salary * employees

voltage = 3.2 # volts
current = 2 # kA / m^2

harvest=700 # $ per kg/hr capacity for Fe

## Tonnage conversions based on stoichiometry
tFe_p_tCl2 = 2 / 3 * 55 / 71 # Roughly 50%, so a 50 tpd iron plant produces 100 tpd cl2
tNaOH_p_tFe = 3*40/55 # Roughly 2.18... 
tFe2O3_p_tFe = 1*160/(2*55.4)
tNaCl_p_tFe = (1/0.6)*3*1000/55.4

## Estimating the capital cost of the stack -----------------------------------

# Specifications on Cu electrowinning blanks from Schlesinger, Mark E.
# Chlorine Specs primarily from Klaucke 2023
Anode_Cost = 1000 # $/m^2 estimated from purchase price of chlor alkali anodes
Cathode_Cost = 250 # $/m^2 estimated from cost of 18 GA SS sheets
Membrane_Cost = 1100 # $/m^2 Hofmann assumes $3000 for a 2.72 m^2 membrane
Ancillary = 3000 # $/m^2 cell area for frames, busbars, flowfields... 
Manufacturing = 2 # cell assembly and profit margins

Stack_Cost = (Anode_Cost + Cathode_Cost + Membrane_Cost + Ancillary) * Manufacturing

# 1980 Membrane + Inert Anode system cost estimated as $1,300 / tpy Cl2 in 2018 dollars
# Stack cost quoted at (E)2250 / kA by Klaucke


## Estimating the capital cost of balance of plant -----------------------------------
# --- Mechanical ---
# Chlorine production per day: 1000 A * 3600s/hr * 24h/day / (2e- * 96485 C/mol) = mol/day
# Chlorine production per year: mol / day * 365 day / year * 70.9 g/mol / (1000 g / kg * 1000 kg / tonne)
Cl2_kA_to_tonne = 1 * 1000 * 3600 * 24 / (2 * 96485) * 365 * 70.9 / (1000 * 1000)

# Values converted from Klaucke from 2021 Euro to 2024 dollars (approx. ~1.1)
Cl2_Treatment = 1388 # $/kgCl2 h–1
NaOH_Storage = 73.7 # $/kgNaOH h–1
Hydrogen_Treatment = 1921 # $/kgH2 h–1
Electrolyte_Circ = 107 # $/kgNaOH h-1 

# --- Electrical ---
Trans_Rect = 562 # $/kA