import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import optimize

def MassBalance(ndot_CO2R_Product_Stream_8 = 41.335978835978835,          #Set Basis for rate of production of reduced product (either CO or C2H4), mol/s (Defaults to 100t C2H4/day)
                Second_Law_Efficiency_Gas_Sep = 0.1,       #Second Law Efficiency of Gas Separation Processes. Assumed to be identical for cathode- and anode-side separation.
                R = 8.314,                                 #Ideal Gas Constant, J/K.mol
                T = 313.15,                                #Temperature, K
                p = 1e5,                                   #Gas Pressure, Pa
                Single_Pass_Conversion_CO2 = 0.125,        #Overall Single Pass Conversion of CO2
                DeltaV = 2.4,                              #Overall Voltage Drop over Electrolyzer, V
                n_CO2R = 12,                               #Number of electrons per mol of CO2R product (12 for C2H4, 2 for CO)
                s_CO2R = 2,                                #Number of CO2 molecules per mol of CO2R product (2 for C2H4, 1 for CO)
                n_H2 = 2,                                  #Number of electrons per mol of H2 produced
                MW_CO2R = 0.028,                           #Molecular Weight, CO2R product, kg/mol
                F = 96485,                                 #Faraday's Constant, C/mol
                Charge_Carbonate_Ion = 2,                  #Charge on Carbonate Ion
                n_OER = 4,                                 #Number of electrons per mol of O2 produced
                Fraction_Charge_Carried_by_CO32 = 1.0,     #Fraction of charge carried by CO32- ions. About 1 for AEM, about 0 for BPM.
                FE_H2_0 = 0.2,                             #Faradaic Efficiency of H2. if Use_Hawks_Model is false, this will just equal the FE of H2. Else, we will take this to be the FE of H2 at very low single pass conversion, when the CO2 concentration is uniform in the reactor and FE_H2 is minimized. 
                i_CO2R_0 = 4000.0,                         #CO2R Partial Current Density, A/m2. Will be used in different ways depending on status of Use_Hawks_Model.
                User_Specify_FE_iCO2R = True,              #If true, use fixed values for electrolyzer performance as specified by user. If false, use Steve Hawk's modela to calcualate how electrolyzer performance changes with single pass conversion. 
                Constant_FE = False,                       #If true, assume FE is constant when using Steve Hawks Model (c.f. 10.1021/acsenergylett.2c01106)
                Constant_iH2 = False,                      #If true, assume iH2 is constant when using Steve Hawks Model (c.f. 10.1021/acsenergylett.2c01106)
                Constant_itot = False):                    #If true, assume itot is constant when using Steve Hawks Model (c.f. 10.1021/acsenergylett.2c01106)

    if [User_Specify_FE_iCO2R, Constant_FE, Constant_iH2, Constant_itot].count(True) != 1:
        raise ValueError('Either User_Specify_FE_iCO2R is True or ONE AND ONLY ONE OF Constant_FE, Constant_iH2 and Constant_itot is True.')
        
    #####################
    #Overall Mass Balance
    #####################
    #The process flow diagram in Figure 1a of the manuscript is sufficiently simple that we can account for the recycle
    #analytically as follows, rather than needing to choose a tear stream and solve iteratively.
    ndot_CO2_Elec_Feed_Stream_2 = ndot_CO2R_Product_Stream_8*s_CO2R/Single_Pass_Conversion_CO2     #Flow rate of CO2 entering electrolyzer, mol/s

    #########################################################
    #Performance of Electrolyzer (Specify Constant FE, iCO2R)
    #########################################################

    #If we don't use Steve Hawks' model (c.f. 10.1021/acsenergylett.2c01106), we will instead use user-specified current density
    if User_Specify_FE_iCO2R is True:
        FE_H2 = FE_H2_0                           #Faradaic Efficiency for H2 Production takes fixed value specified by user.
        i_CO2R = i_CO2R_0                         #CO2R Partial Current Density of electrolyzer takes fixed value specified by user. A/m2
        i_H2 = i_CO2R * FE_H2/(1-FE_H2)            #H2 Partial Current Density, A/m2
        i_tot = i_CO2R + i_H2                     #Total Current Density, A/m2
        A = ndot_CO2_Elec_Feed_Stream_2*Single_Pass_Conversion_CO2*n_CO2R*F/(s_CO2R*i_CO2R)       #Electrolyzer Area, m2
        eta = None
       
    ############################################
    #Performance of Electrolyzer (Hawks, Fix FE)
    ############################################

    #We use Steve Hawks model (c.f. 10.1021/acsenergylett.2c01106), assuming constant FE.
    elif (not User_Specify_FE_iCO2R and Constant_FE) is True:

        eta = Single_Pass_Conversion_CO2 * (1 + n_CO2R*Fraction_Charge_Carried_by_CO32/((1-FE_H2_0)*s_CO2R*Charge_Carbonate_Ion))    #Single Pass CO2 Consumption
        i_CO2R = -i_CO2R_0 * eta/np.log(1-eta)                  #CO2R Partial Current Density, A/m2
        i_H2 = i_CO2R * FE_H2_0/(1-FE_H2_0)                     #H2 Partial Current Density, A/m2
        i_tot = i_CO2R + i_H2                                   #Total Current Density, A/m2
        FE_H2 = i_H2/i_tot 
        A = ndot_CO2_Elec_Feed_Stream_2*Single_Pass_Conversion_CO2*n_CO2R*F/(s_CO2R*i_CO2R)     #Electrolyzer Area, m2

    #############################################
    #Performance of Electrolyzer (Hawks, Fix iH2)
    #############################################

    #We use Steve Hawks model (c.f. 10.1021/acsenergylett.2c01106) to estimate i_CO2R, but assume iH2 is independent of gas CO2 concentration.
    elif (not User_Specify_FE_iCO2R and Constant_iH2) is True:

        eta = Single_Pass_Conversion_CO2 * (1 + n_CO2R*Fraction_Charge_Carried_by_CO32/((1-FE_H2_0)*s_CO2R*Charge_Carbonate_Ion))    #Single Pass CO2 Consumption
        i_CO2R = -i_CO2R_0*eta/np.log(1-eta)                    #CO2R Partial Current Density, A/m2
        i_H2 = i_CO2R_0 * FE_H2_0/(1-FE_H2_0)                   #H2 Partial Current Density, A/m2
        i_tot = i_CO2R + i_H2                                   #Total Current Density, A/m2
        A = ndot_CO2_Elec_Feed_Stream_2*Single_Pass_Conversion_CO2*n_CO2R*F/(s_CO2R*i_CO2R)     #Electrolyzer Area, m2
        FE_H2 = i_H2/i_tot     
        
    ##############################################
    #Performance of Electrolyzer (Hawks, Fix itot)
    ##############################################

    #We use Steve Hawks model (c.f. 10.1021/acsenergylett.2c01106) to estimate i_CO2R, but assume itot is independent of gas CO2 concentration.
    elif (not User_Specify_FE_iCO2R and Constant_itot) is True:

        eta = Single_Pass_Conversion_CO2 * (1 + n_CO2R*Fraction_Charge_Carried_by_CO32/((1-FE_H2_0)*s_CO2R*Charge_Carbonate_Ion))                          #Single Pass CO2 Consumption
        i_CO2R = -i_CO2R_0 * eta/np.log(1-eta)                  #CO2R Partial Current Density, A/m2
        i_tot = i_CO2R_0 / (1-FE_H2_0)
        i_H2 = i_tot - i_CO2R
        A = ndot_CO2_Elec_Feed_Stream_2*Single_Pass_Conversion_CO2*n_CO2R*F/(s_CO2R*i_CO2R)     #Electrolyzer Area, m2
        FE_H2 = i_H2/i_tot    
     
   
    ####################################
    #FLOW INTO ANODE-SIDE GAS SEPARATION
    ####################################

    ndot_O2_Anode_Feed_Stream_4 = i_tot*A/(n_OER*F)
    ndot_CO2_Anode_Feed_Stream_4 = i_tot*A/(Charge_Carbonate_Ion*F)*Fraction_Charge_Carried_by_CO32

    ######################################
    #FLOW INTO CATHODE-SIDE GAS SEPARATION
    ######################################

    ndot_H2_Cathode_Feed_Stream_3 = i_H2*A/(n_H2*F)
    ndot_CO2_Cathode_Feed_Stream_3 = ndot_CO2_Elec_Feed_Stream_2 \
                                     - ndot_CO2_Elec_Feed_Stream_2*Single_Pass_Conversion_CO2 \
                                     - ndot_CO2_Anode_Feed_Stream_4
    ndot_CO2R_Cathode_Feed_Stream_3 = ndot_CO2_Elec_Feed_Stream_2*Single_Pass_Conversion_CO2/s_CO2R
    
    #############################
    #ELECTROLYZER ENERGY COST (W)
    #############################

    E_Electrolyzer = DeltaV*i_tot*A        

    #################################
    #ANODE SEPARATION ENERGY COST (W)
    #################################

    ndot_tot_Anode = ndot_CO2_Anode_Feed_Stream_4 + ndot_O2_Anode_Feed_Stream_4
    y_O2_Anode = ndot_O2_Anode_Feed_Stream_4 / ndot_tot_Anode
    y_CO2_Anode = ndot_CO2_Anode_Feed_Stream_4 / ndot_tot_Anode
    E_Anode_Sep = -ndot_tot_Anode*R*T*(y_O2_Anode*np.log(y_O2_Anode) \
                + y_CO2_Anode*np.log(y_CO2_Anode))/Second_Law_Efficiency_Gas_Sep

    ###################################
    #CATHODE SEPARATION ENERGY COST (W)
    ###################################

    ndot_tot_Cathode = ndot_CO2_Cathode_Feed_Stream_3 + \
                        ndot_H2_Cathode_Feed_Stream_3 + ndot_CO2R_Cathode_Feed_Stream_3

    y_H2_Cathode = ndot_H2_Cathode_Feed_Stream_3 / ndot_tot_Cathode
    y_CO2_Cathode = ndot_CO2_Cathode_Feed_Stream_3 / ndot_tot_Cathode
    y_CO2R_Cathode = ndot_CO2R_Cathode_Feed_Stream_3 / ndot_tot_Cathode

    E_Cathode_Sep = -ndot_tot_Cathode*R*T\
                    *((y_CO2R_Cathode+y_H2_Cathode)*np.log(y_H2_Cathode+y_CO2R_Cathode)\
                     +y_CO2_Cathode*np.log(y_CO2_Cathode))/Second_Law_Efficiency_Gas_Sep       #J

    sol = {
        "E_Electrolyzer": E_Electrolyzer,
        "E_Anode_Sep": E_Anode_Sep,
        "E_Cathode_Sep": E_Cathode_Sep,
        "E_total" : E_Electrolyzer + E_Anode_Sep + E_Cathode_Sep,
        "ndot_CO2_Anode_Feed_Stream_4": ndot_CO2_Anode_Feed_Stream_4,
        "ndot_O2_Anode_Feed_Stream_4": ndot_O2_Anode_Feed_Stream_4,
        "ndot_CO2_Cathode_Feed_Stream_3": ndot_CO2_Cathode_Feed_Stream_3,
        "ndot_H2_Cathode_Feed_Stream_3": ndot_H2_Cathode_Feed_Stream_3,
        "ndot_CO2R_Cathode_Feed_Stream_3": ndot_CO2R_Cathode_Feed_Stream_3,
        "ndot_CO2R_Product_Stream_8": ndot_CO2R_Product_Stream_8,
        "ndot_CO2_Elec_Feed_Stream_2": ndot_CO2_Elec_Feed_Stream_2,        
        "A": A,
        "i_H2": i_H2,    
        "i_CO2R": i_CO2R,
        "eta": eta,
        "i_tot": i_tot,
        "FE_H2": FE_H2
    }
    
    return sol


def ElectrolyzerCosting(
    Convertunit = 24*3600/1000,                #Convert kg/s to t/d
    MW_O2 = 0.032,                             #Molecular weight of O2, kg/mol
    MW_CO2 = 0.044,                            #Molecular weight of CO2, kg/mol
    MW_H2 = 0.002,                             #Molecular weight of H2, kg/mol
    Plantlife = 20.0,                          #Plant Life, years
    operatingdays = 328.5,                     #Operating days, days/year
    ElectricityPrice = 0.02,                   #Electricity Price, $/kWh
    PEMCost = 1250.0,                          #Electrolyzer Cost, $/kW
    CAPEXGasSepReference = 1989043.0,          #PSA CAPEX Reference, for PSA capacity of 1000 m3/hr, $
    PSACapacityReference = 1000.0,             #PSA Capacity Reference, m3/hr
    PSAOPEX = 0.25,                            #PSA OPEX, kwh/m3
    PSAScaleFactor = 0.7,                      #PSA Scale Factor
    ndot_CO2R_Product_Stream_8 = 41.335978835978835,          #Set Basis for rate of production of reduced product (either CO or C2H4), mol/s (Defaults to 100t C2H4/day)
    Second_Law_Efficiency_Gas_Sep = 0.1,       #Second Law Efficiency of Gas Separation Processes. Assumed to be identical for cathode- and anode-side separation.
    R = 8.314,                                 #Ideal Gas Constant, J/K.mol
    T = 313.15,                                #Temperature, K
    p = 1e5,                                   #Gas Pressure, Pa
    Single_Pass_Conversion_CO2 = 0.125,        #Overall Single Pass Conversion of CO2
    DeltaV = 2.4,                              #Overall Voltage Drop over Electrolyzer, V
    n_CO2R = 12,                               #Number of electrons per mol of CO2R product (12 for C2H4, 2 for CO)
    s_CO2R = 2,                                #Number of CO2 molecules per mol of CO2R product (2 for C2H4, 1 for CO)
    n_H2 = 2,                                  #Number of electrons per mol of H2 produced
    MW_CO2R = 0.028,                           #Molecular Weight, CO2R product, kg/mol
    F = 96485,                                 #Faraday's Constant, C/mol
    Charge_Carbonate_Ion = 2,                  #Charge on Carbonate Ion
    n_OER = 4,                                 #Number of electrons per mol of O2 produced
    Fraction_Charge_Carried_by_CO32 = 1.0,     #Fraction of charge carried by CO32- ions. About 1 for AEM, about 0 for BPM.
    FE_H2_0 = 0.2,                             #Faradaic Efficiency of H2. if Use_Hawks_Model is false, this will just equal the FE of H2. Else, we will take this to be the FE of H2 at very low single pass conversion, when the CO2 concentration is uniform in the reactor and FE_H2 is minimized. 
    i_CO2R_0 = 4000.0,                         #CO2R Partial Current Density, A/m2. Will be used in different ways depending on status of Use_Hawks_Model.
    User_Specify_FE_iCO2R = True,              #If true, use fixed values for electrolyzer performance as specified by user. If false, use Steve Hawk's modela to calcualate how electrolyzer performance changes with single pass conversion. 
    Constant_FE = False,                       #If true, assume FE is constant when using Steve Hawks Model
    Constant_iH2 = False,                      #If true, assume iH2 is constant when using Steve Hawks Model
    Constant_itot = False,                     #If true, assume itot is constant when using Steve Hawks Model
    Jouny_PSA_OPEX = True,                     #If true, use Jouny's PSA OPEX correlation (equivalent to 2nd law efficiency of ~7%). If false, use second law efficiency directly.
    Separate_Plant_Lives = False,              #If true, calculate electrolyzer and PSA plant life seaprately
    Electrolyzer_Life = 20,                    #Used for electrolyzer lifetime if Separate_Plant_Lives = True
    Fixed_sep_cost_per_tCO2 = False,           #If True, use fixed cost of separation per tCO2 (e.g. $40/tCO2)
    Sep_cost_per_tCO2 = 25,                    #Cost of CO2 separation if Fixed_sep_cost_per_tCO2 = True
    Wilcox_sep_cost_per_tCO2 = False):

                    
    #Solve Mass Balance
    MB = MassBalance(ndot_CO2R_Product_Stream_8 = ndot_CO2R_Product_Stream_8,
        Second_Law_Efficiency_Gas_Sep = Second_Law_Efficiency_Gas_Sep,
        R = R,
        T = T,
        p = p,
        Single_Pass_Conversion_CO2 = Single_Pass_Conversion_CO2,
        DeltaV = DeltaV,
        n_CO2R = n_CO2R,
        s_CO2R = s_CO2R,
        n_H2 = n_H2,
        MW_CO2R = MW_CO2R,
        F = F,
        Charge_Carbonate_Ion = Charge_Carbonate_Ion,
        n_OER = n_OER,
        Fraction_Charge_Carried_by_CO32 = Fraction_Charge_Carried_by_CO32,
        FE_H2_0 = FE_H2_0,
        i_CO2R_0 = i_CO2R_0,
        User_Specify_FE_iCO2R = User_Specify_FE_iCO2R,
        Constant_FE = Constant_FE,
        Constant_iH2 = Constant_iH2,
        Constant_itot = Constant_itot)
    
    #Calculate Gas Density (kg/m3)
    CO2RDensity = p/(R*T)*MW_CO2R                #CO2R Product Density, kg/m3
    H2Density = p/(R*T)*MW_H2                    #Hydrogen Density, kg/m3
    CO2Density = p/(R*T)*MW_CO2                   #CO2 Density, kg/m3
    O2Density = p/(R*T)*MW_O2                    #O2 Density, kg/m3

   # CO2RDensity = 1.18 #kg/m3
   # H2Density = 0.082#kg/m3
   # CO2Density = 1.98 #kg/m3
   # O2Density = 1.429 #kg/m3

    #Economic Calculations
    FE_H2 = MB['FE_H2']
    A = MB['A']
    Total_Energy = MB['E_total']
    i_CO2R = MB['i_CO2R']
    i_tot = MB['i_tot']
    E_Electrolyzer = MB['E_Electrolyzer']
    E_Anode_Sep = MB['E_Anode_Sep']
    E_Cathode_Sep = MB['E_Cathode_Sep']
    Anode_gas_vol = MB['ndot_O2_Anode_Feed_Stream_4']*MW_O2*Convertunit*1000/O2Density/24+MB['ndot_CO2_Anode_Feed_Stream_4']*MW_CO2*Convertunit*1000/CO2Density/24    
    Cathode_gas_vol = MB['ndot_H2_Cathode_Feed_Stream_3']*MW_H2*Convertunit*1000/H2Density/24+MB['ndot_CO2_Cathode_Feed_Stream_3']*MW_CO2*Convertunit*1000/CO2Density/24+MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit*1000/CO2RDensity/24    
    if Separate_Plant_Lives:
        Electrolyzercapexperproduct = PEMCost*3.6*A*1.3/Electrolyzer_Life/(MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)/operatingdays                       #$/t C2H4
    else:
        Electrolyzercapexperproduct = PEMCost*3.6*A*1.3/Plantlife/(MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)/operatingdays                                                                                                                      #$/t C2H4 
    Electrolyzeropexperproduct = ElectricityPrice*(MB['E_Electrolyzer']/1000*24)/(MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)
    PSAcapexperproduct = (CAPEXGasSepReference*(Anode_gas_vol/PSACapacityReference)**PSAScaleFactor+CAPEXGasSepReference*(Cathode_gas_vol/PSACapacityReference)**PSAScaleFactor)/Plantlife/(MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)/operatingdays 
    if Jouny_PSA_OPEX:
        PSAopexperproduct = (ElectricityPrice*PSAOPEX*(Anode_gas_vol + Cathode_gas_vol)*24)/(MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)
    else:
        PSAopexperproduct = ElectricityPrice*((E_Anode_Sep+E_Cathode_Sep)/1000*24)/(MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)
        
    Separation_Cost = PSAcapexperproduct + PSAopexperproduct
    if Fixed_sep_cost_per_tCO2:
        Separation_Cost = Sep_cost_per_tCO2*((MB['ndot_CO2_Cathode_Feed_Stream_3'] + MB['ndot_CO2_Anode_Feed_Stream_4'])*MW_CO2*Convertunit) / (MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)
       
    if Wilcox_sep_cost_per_tCO2:
        α_Anode = 0.95      #CO2 capture fraction
        α_Cathode = 0.95    #CO2 capture fraction
        β_Anode =  MB['ndot_CO2_Anode_Feed_Stream_4'] / ( MB['ndot_CO2_Anode_Feed_Stream_4'] +  MB['ndot_O2_Anode_Feed_Stream_4']) * 100
        β_Cathode = MB['ndot_CO2_Cathode_Feed_Stream_3'] / ( MB['ndot_CO2_Cathode_Feed_Stream_3'] +  MB['ndot_H2_Cathode_Feed_Stream_3'] + MB['ndot_CO2R_Cathode_Feed_Stream_3']) * 100
        Sep_cost_per_tCO2_Anode   = 10**(2.18 - 0.426*np.log10(α_Anode)   - 0.391*np.log10(β_Anode)   - 0.028*np.log10(1000))
        Sep_cost_per_tCO2_Cathode = 10**(2.18 - 0.426*np.log10(α_Cathode) - 0.391*np.log10(β_Cathode) - 0.028*np.log10(1000))
        Separation_Cost_Cathode = Sep_cost_per_tCO2_Cathode*((MB['ndot_CO2_Cathode_Feed_Stream_3'])*MW_CO2*Convertunit) / (MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)
        Separation_Cost_Anode = Sep_cost_per_tCO2_Anode*((MB['ndot_CO2_Anode_Feed_Stream_4'])*MW_CO2*Convertunit) / (MB['ndot_CO2R_Product_Stream_8']*MW_CO2R*Convertunit)
        Separation_Cost = Separation_Cost_Anode + Separation_Cost_Cathode
    Electrolyzer_Cost = Electrolyzercapexperproduct + Electrolyzeropexperproduct
    Total_Cost = Separation_Cost + Electrolyzer_Cost      
        
    sol = {
        "Anode_gas_vol": Anode_gas_vol,
        "Cathode_gas_vol": Cathode_gas_vol,
        "Electrolyzercapexperproduct": Electrolyzercapexperproduct, 
        "Electrolyzeropexperproduct": Electrolyzeropexperproduct,
        "PSAcapexperproduct": PSAcapexperproduct,
        "PSAopexperproduct": PSAopexperproduct,
        "Separation_Cost":Separation_Cost,
        "Electrolyzer_Cost": Electrolyzer_Cost,
        "Total_Cost": Total_Cost
    }
    
    sol.update(MB)
    
    return sol
    