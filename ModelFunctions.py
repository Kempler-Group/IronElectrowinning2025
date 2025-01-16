'''
This file contains functions required to reproduce figures presented in Konovalova et al 2024
'''

import pandas as pd

## Reference list
#1. Klaucke, https://pubs.acs.org/doi/full/10.1021/acs.iecr.2c04188
#2. Hofmann, https://aiche.onlinelibrary.wiley.com/doi/epdf/10.1002/aic.17480
#3. Lahrsen, https://www.mdpi.com/2227-9717/10/4/761
#4. Stinn, https://iopscience.iop.org/article/10.1149/2.F06202IF/pdf

def stack_cost(CURRENT = 2, alkaline=False):
    
    from ModelAssumptions import Membrane_Cost, Anode_Cost, Cathode_Cost, Ancillary, Manufacturing

    # Current in kA / m^2
    # All costs defined above in $/m^2
    # Returns cost per kA
    
    # add a flag for removing membrane costs
    
    if alkaline:
        MEMBRANE = 0
    else:
        MEMBRANE = Membrane_Cost
    
    return (Anode_Cost + Cathode_Cost + MEMBRANE + Ancillary) * Manufacturing / CURRENT


def BOP_cost(Fe_FE = 0.9, NaOH_FE = 0.95, Cl2_FE = 0.98, Fe_Harvest = 700, printout=False, reportout=False):
    
    from ModelAssumptions import Cl2_Treatment, NaOH_Storage, Hydrogen_Treatment, Electrolyte_Circ, Trans_Rect 
    
    # Cost of harvesting iron is not proven and should be swept as a parameter ($ per kg/hr capacity)
    # We have assumed the cost to lie somewhere in between the cost of processing Cl2 and the cost of processing NaOH
    
    # Supporting figure - cost of iron harvesting
    
    # Function returns the cost per system kA
    basis_current = 1000 # Amperes
    
    F = 96485 # C / mol
    
    H2_FE = 1 - Fe_FE
    Rate_H2 = basis_current * H2_FE * 3600 / (2 * F) * 2 / 1000 # kgH2 / hr
    Rate_Cl2 = basis_current * Cl2_FE * 3600 / (2 * F) * 70.9 / 1000 # kgCl2 / hr
    Rate_Fe = basis_current * Fe_FE * 3600 / (3 * F) * 55.85 / 1000 # kgFe / hr
    Rate_NaOH = basis_current * NaOH_FE * 3600 / (F) * 39.9 / 1000 # kgNaOH / hr
    
    H2_Cost = Hydrogen_Treatment * Rate_H2
    Cl2_Cost = Cl2_Treatment * Rate_Cl2
    NaOH_Cost = (NaOH_Storage + Electrolyte_Circ) * Rate_NaOH
    Fe_Cost = Fe_Harvest * Rate_Fe
    
    if printout:
        print("Balance of plant costs, $/kA system")
        print("Hydrogen: ", H2_Cost)
        print("Chlorine: ", Cl2_Cost)
        print("NaOH: ", NaOH_Cost)
        print("Iron: ", Fe_Cost)
        
    if reportout:
        
        costs = {'Fe': Fe_Cost,
                 'NaOH': NaOH_Cost,
                 'Cl2': Cl2_Cost,
                 'H2': H2_Cost,
                 'Electrical':Trans_Rect
        }
        
        return costs
    
    return Fe_Cost + NaOH_Cost + Cl2_Cost + H2_Cost + Trans_Rect
    

def alkaline_BOP_cost(Fe_FE = 0.9, NaOH_FE = 0.95, Fe_Harvest = 700, printout=False, reportout=False):

    from ModelAssumptions import Hydrogen_Treatment, Electrolyte_Circ, Trans_Rect

    # For comparison with the primary chlor iron model
    # Function returns the cost per system kA
    basis_current = 1000 # Amperes
    
    F = 96485 # C / mol
    
    H2_FE = 1 - Fe_FE
    Rate_H2 = basis_current * H2_FE * 3600 / (2 * F) * 2 / 1000 # kgH2 / hr
    Rate_Fe = basis_current * Fe_FE * 3600 / (3 * F) * 55.85 / 1000 # kgFe / hr
    
    # Assume similar recirculation rates to chlor-iron stack
    Rate_NaOH = basis_current * 3600 / (F) * 39.9 / 1000 # kgNaOH / hr
    
    #Removed Cl2 costs, still need H2 treatment
    H2_Cost = Hydrogen_Treatment * Rate_H2
    #Assumed cost for recirculation of NaOH but no storage necessary
    #Not included are the additional costs necessary for NaOH treatment
    NaOH_Cost = (Electrolyte_Circ) * Rate_NaOH
    Fe_Cost = Fe_Harvest * Rate_Fe
    
    if printout:
        print("Balance of plant costs, $/kA system")
        print("Hydrogen: ", H2_Cost)
        print("NaOH: ", NaOH_Cost)
        print("Iron: ", Fe_Cost)
        
    if reportout:
        
        costs = {'Fe': Fe_Cost,
                 'NaOH': NaOH_Cost,
                 'H2': H2_Cost,
                 'Electrical':Trans_Rect
        }
        
        return costs
    
    return Fe_Cost + NaOH_Cost +  H2_Cost + Trans_Rect


def Fe_kA_to_tpy(Fe_FE=0.9, basis_current = 1, n_e = 3):
    
    # Calculate conversion from cell kA to tpy Fe production
    # tonne per year Fe production for the given basis current (in kA)
    
    F = 96485 # C/mol
    
    # Convert kA to A and then multiply by number of electrons, faraday's constant, g/kg, iron molar mass
    Rate_Fe = basis_current * 1000 * Fe_FE * 3600 / (n_e * F) * 55.85 / 1000 # kgFe / hr
    
    # Return hourly rate (kg/hr) of Fe times 8760 (hours/year) / 1000 (kg/tonne)
    
    return Rate_Fe * 8760 / 1000


def pol_curve(J, ASR = 5, E0 = 2.2):
    
    # ASR in ohm-cm^2, E0 in volts
    # this function takes a current density as an input (kA / m^2)
    # then converts to A / cm^2
    # the function returns the cell voltage as an output
    J = J / 10 # convert to A / cm^2 
    
    E = E0 + J * ASR
    
    return E

def efficiency_fe(volts, Fe_FE=1):
    ## calculates efficiency in Wh/g, kWh/kg, MWh/tonne
    eff = volts*3*96485/(3600*55.4)/Fe_FE
    return eff

import ModelAssumptions as BASE

def NPV_calc(Iron_Price, Iron_Prod=BASE.prod_tpd, Cell_Voltage = BASE.voltage, 
             Electric_Price=BASE.EP, CellCurrent=BASE.current, Fe_Selectivity = 0.9,
             CAPACITY=0.98, NaOHCAPACITY=0.95, Cl2CAPACITY=0.95, 
             iron_ore_price=120, Harvest=BASE.harvest,
             brine_price=1, naoh_price=0, cl2_price=BASE.Cl2P, replace_rate=7,
             chloriron_cell=True, alkaline_cell=False, PolCurve=False, ASR=5):
    
    # Default model assumes no NaOH sales, 7 year lifetime for the cells in the stack
    # Units of current are kA/m^2
    # ASR is area specific resistance in ohms-cm^2
    
    #Timeline
    model_years = BASE.life #years
    
    #Scheduled cash flow
    CashFlow = pd.DataFrame({'Years':[],'Stack Life':[],'Inflation Year':[],
                            'Stack Eff (MWh/t)':[], 'Sales':[],'Replacement Costs':[],'Operating Costs':[],
                            'Net Cash Flow':[],'Discounted Flow':[]})
    
    if chloriron_cell:
        if alkaline_cell:
            raise Exception('Select a single cell class')
            
        StackCost_kA = stack_cost(CURRENT = CellCurrent, alkaline=False) # Operating current in kA / m^2
        BoPCost_kA = BOP_cost(Fe_FE = 0.9, NaOH_FE = 0.95, Cl2_FE = 0.98, Fe_Harvest = Harvest) # Harvest cost 
    
    elif alkaline_cell:
        StackCost_kA = stack_cost(CURRENT = CellCurrent, alkaline=True) # Operating current in kA / m^2
        BoPCost_kA = alkaline_BOP_cost(Fe_FE = 0.9, NaOH_FE = 0.95, Fe_Harvest = Harvest) # Harvest cost
        Cl2CAPACITY=0
    
    else:
        raise Exception('Select a cell type for estimating stack and BoP')
    
    #Keep voltage drift = 0 for now
    #voltage_drift = 0
    
    #Calculate Stack Eff
    #If pol curve flag is on, calculate the cell voltage based on the specified current (kA m^-2)
    
    if chloriron_cell:
        Initial_Voltage = pol_curve(CellCurrent, ASR=ASR, E0=2.2) if PolCurve else Cell_Voltage #else user specified voltage
    
    elif alkaline_cell:
        Initial_Voltage = pol_curve(CellCurrent, ASR=ASR, E0=1.4) if PolCurve else Cell_Voltage #else user specified voltage
    
    else:
        raise Exception('Select a cell type for estimating cell voltage')
    
    stackeff = efficiency_fe(Initial_Voltage, Fe_FE=Fe_Selectivity) # MWh/tonne
    
    # Convert to tpy and multiply by nameplate plant capacity
    StackCost = StackCost_kA / Fe_kA_to_tpy(Fe_FE=Fe_Selectivity) * Iron_Prod * 365  
    BoPCost = BoPCost_kA / Fe_kA_to_tpy(Fe_FE=Fe_Selectivity) * Iron_Prod * 365
    SystemCost = StackCost + BoPCost # this is in dollars
    
    # For troubleshooting
    # print("\nStack Cost = {0:.2e}".format(StackCost))
    # print("BoP Cost = {0:0.2e}".format(BoPCost))
    # print("System Cost = {0:0.2e}".format(SystemCost))
    
    ## Initialize variables
    year = [1]
    stacklife = [0]
    inflationyear = [1]
    sales = [0]
    replacement_costs = [0]
    operating_costs = [0]
    net_cash = [-SystemCost]
    discounted_cash = [net_cash[-1]/(1+BASE.IRR)**year[-1]]
    
    for x in range(2,model_years+2):
        year.append(x)
        inflationyear.append((1+BASE.inflation)**x)

        # Calculate Iron, Chlorine, NaOH sales
        iron_sale = Iron_Price*Iron_Prod*365*CAPACITY
        cl2_sale = cl2_price*Iron_Prod/BASE.tFe_p_tCl2*365*CAPACITY*Cl2CAPACITY
        naoh_sale = naoh_price*Iron_Prod*BASE.tNaOH_p_tFe*365*CAPACITY*NaOHCAPACITY
        
        # Track sales
        sales.append(iron_sale+cl2_sale+naoh_sale) # [$/kg] * [kg/h] * [h/year] * percent
        sales[-1] = sales[-1]*inflationyear[-1]

        # Track the stack life and repurchase stack if year > stack life
        if stacklife[-1] < replace_rate and year[-1] != model_years: # Don't replace stack in the last year of life
            stacklife.append(stacklife[-1]+1)
        else:
            stacklife.append(1)

        #Maintenance cost are 0.5% system cost per year + stack replacement
        replacement_costs.append(StackCost*(stacklife[-1]==replace_rate) + 0.005*SystemCost)
        replacement_costs[-1] = replacement_costs[-1]*inflationyear[-1]

        #Operating costs are assumed to be dominated by electricity prices and ore prices
        #Also included are water expenses
        ann_prod = Iron_Prod*365*CAPACITY # iron tpy
        e_expense = ann_prod*stackeff*Electric_Price # $/MWh * MWh/tonne * tpy
        w_expense = brine_price*ann_prod*BASE.tNaCl_p_tFe
        fe2o3_expense = iron_ore_price*ann_prod*BASE.tFe2O3_p_tFe
        
        if Iron_Prod == 100:
            annual_labor = BASE.labor
        else:
            raise Exception("Warning, check plant scale and labor costs")
        
        expenses = e_expense+w_expense+fe2o3_expense+annual_labor
        #[$/kWh] * [kWh/kg] * [kg/h] * [h/year] * percent
        operating_costs.append(expenses)
        operating_costs[-1] = operating_costs[-1]*inflationyear[-1]

        #Net cash flow = Sales - replacement - operating costs
        net_cash.append(sales[-1] - operating_costs[-1] - replacement_costs[-1])

        #Discounted cash flow = net cash flow / (1+IRR)^year
        discounted_cash.append(net_cash[-1]/(1+BASE.IRR)**year[-1])
        
    CashFlow['Years'] = year
    CashFlow['Stack Life'] = stacklife
    CashFlow['Inflation Year'] = inflationyear
    CashFlow['Stack Eff (MWh/t)'] = stackeff
    CashFlow['Sales'] = sales
    CashFlow['Replacement Costs'] = replacement_costs
    CashFlow['Operating Costs'] = operating_costs
    CashFlow['Net Cash Flow'] = net_cash
    CashFlow['Discounted Flow'] = discounted_cash
    
    return sum(CashFlow['Discounted Flow']), CashFlow

def LCOFe(Iron_Prod=BASE.prod_tpd, Cell_Voltage = BASE.voltage, 
         Electric_Price=BASE.EP, CellCurrent=BASE.current, Fe_Selectivity = 0.9,
         CAPACITY=0.98, NaOHCAPACITY=0.95, Cl2CAPACITY=0.95, 
         iron_ore_price=120, Harvest=BASE.harvest,
         brine_price=1, naoh_price=0, cl2_price=BASE.Cl2P, replace_rate=7,
         chloriron_cell=True, alkaline_cell=False, PolCurve=False, ASR=5):
    
    low_price = 1
    high_price = 1000
    
    low_val, low_table = NPV_calc(low_price, Iron_Prod, Cell_Voltage, 
                         Electric_Price, CellCurrent, Fe_Selectivity,
                         CAPACITY, NaOHCAPACITY, Cl2CAPACITY, 
                         iron_ore_price, Harvest,
                         brine_price, naoh_price, cl2_price, replace_rate,chloriron_cell,alkaline_cell,
                         PolCurve, ASR)
    
    high_val, high_table = NPV_calc(high_price, Iron_Prod, Cell_Voltage, 
                             Electric_Price, CellCurrent, Fe_Selectivity,
                             CAPACITY, NaOHCAPACITY, Cl2CAPACITY, 
                             iron_ore_price, Harvest,
                             brine_price, naoh_price, cl2_price, replace_rate,chloriron_cell,alkaline_cell,
                             PolCurve, ASR)
    
    # Interpolate to solve for exact LCOX
    slope = (high_price - low_price)/(high_val - low_val)

    return high_price - high_val*slope