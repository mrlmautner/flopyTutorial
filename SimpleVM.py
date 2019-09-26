# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 16:29:31 2019

@author: MM
"""

import flopy
import numpy as np
import calendar

def initializeFM(modelname, xll, yll, xur, yur, cellsize, STRT_YEAR, END_YEAR, ACTIVE, GEO, DEM, ZoneParams, IH='DEM', THICK1=50, THICK2=350, nlay=2):
    
    # modelname to set the file root 
    mf = flopy.modflow.Modflow(modelname, exe_name=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe')
    
    # Model domain and grid definition
    ztop = DEM # Model top elevation (Digital Elevation Model)
    L1botm = ztop - THICK1 # Layer 1 bottom elevation
    L2botm = L1botm - THICK2 # Layer 2 bottom elevation
    botm = [L1botm,L2botm] # Model bottom elevation
    ncol = int((xur-xll)/cellsize) # Number of rows
    nrow = int((yur-yll)/cellsize) # Number of columns
    delr = cellsize # Row height
    delc = cellsize # Column height
    
    # Time discretization
    nper = (END_YEAR - STRT_YEAR)*12 # Number of stress periods
    nstp = []
    for y in range(STRT_YEAR,END_YEAR):
        for m in range(1,13):
            nstp.append(calendar.monthrange(y,m)[1])
    nstp = np.array(nstp) # Assign the number of time steps in each stress period
    steady = np.zeros((nper),dtype=bool) # Assign stress period as steady state or transient
    
    dis = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, delr=delr, delc=delc, top=ztop, botm=botm, perlen=nstp, nstp=nstp, steady=steady, start_datetime='01/01/1984')
        
    # Model Boundaries & initial conditions
    # Active areas: Multiply active array of 1s and 0s by array of 1s
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    for a, activearray in enumerate(ACTIVE):
        ibound[a,:,:] = ibound[a,:,:]*activearray
    
    # Variables for the BAS package
    if IH is 'DEM':
        strt = DEM # Assign model top as initial head
    else:
        strt = IH
    
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=3)
    
    # Layer properties
    # Add LPF package to the MODFLOW model
    HK = np.zeros((nlay,nrow,ncol)) # Hydraulic conductivity
    VKA = np.zeros((nlay,nrow,ncol)) # Vertical anisotropy (H:V) of hydraulic conductivity
    SS = np.zeros((nlay,nrow,ncol)) # Specific storage
    SY = np.zeros((nlay,nrow,ncol)) # Specific yield
    
    # Loop through the layers and formations for each layer to apply the geologic parameters to each array
    for l in range(nlay):
        
        for f in ZoneParams['HK'][0][l]:
            HK[l,:,:] += (GEO[l] == f) * ZoneParams['HK'][1][f]
        
        # Vertical anisotropy (H:V) of hydraulic conductivity
        for f in ZoneParams['VKA'][0][l]:
            VKA[l,:,:] += (GEO[l] == f) * ZoneParams['VKA'][1][f]
        
        # Specific storage
        for f in ZoneParams['SS'][0][l]:
            SS[l,:,:] += (GEO[l] == f) * ZoneParams['SS'][1][f]
        
        # Specific yield
        for f in ZoneParams['SY'][0][l]:
            SY[l,:,:] += (GEO[l] == f) * ZoneParams['SY'][1][f]
    
    layvka = [1]*nlay # Indicates that VKA represents the ratio of H:V hydraulic conductivity
    
    lpf = flopy.modflow.mflpf.ModflowLpf(mf, layvka=layvka, hk=HK, vka=VKA, ss=SS, sy=SY)
    
    return ncol, nrow, mf, dis, bas, lpf

def addRecharge(LU_arrays, PRECIP, RCH_mult, yrzero=0, S_YR=0, E_YR=0, RCH_Dict=0):
    # LU_arrays: dictionary with an entry for each land use type which contains gridded percent amounts for each land use type
    # PRECIP: dictionary with an entry for each stress period which contains gridded precipitation
    # RCH_Dict: existing dictionary holding recharge data or 0 if the dictionary must be initialized
    
    # Initialize dictionary: if there is no exisiting dictionary, create dictionary with no entries
    if RCH_Dict == 0:
        RCH_Dict = {}
    
    # If the recharge is for the first time step, apply only to the first time step
    if S_YR == 0:
        for l, landuse in enumerate(['URBAN','NATURAL','WATER']):
            # If there is not already an entry for the selected stress period, create a new array
            try:
                RCH_Dict[0] += PRECIP[0]*LU_arrays[landuse]*RCH_mult[l]
            except:
                RCH_Dict[0] = PRECIP[0]*LU_arrays[landuse]*RCH_mult[l]
    
    # Loop through all stress periods between S_YR and E_YR
    else:
        for year in range(int(S_YR),int(E_YR)):
            for month in range(0,12):
                
                per = (year - yrzero)*12 + month
                
                # Apply recharge amounts for each land use type                
                for l, landuse in enumerate(['URBAN','NATURAL','WATER']):                    
                    # If there is not already an entry for the selected stress period, create a new array
                    try:
                        RCH_Dict[per] += PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]
                    except:
                        RCH_Dict[per] = PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]

    return RCH_Dict

def addNewWells(New_WEL, LYR, xll, yur, cellsize, WEL_Dict=0, WEL_mult=1, coordType='rc'):
    # New_WEL is an np array of the following format: X (or C), Y (or R), Start Period, End Period, Flow (m3/t)
    # WEL_PAR is a scalar multiplier to be applied to all wells in the data set New_WEL
    # WEL_Dict is a dictionary that contains dictionary for each stress period, each dictionary contains an entry for each well with the layer, row, column, and pumping rate
    # coordType is a marker that should be either 'xy' or 'rc' depending on the coordinate definition in the New_WEL array
    
    # Initialize dictionary    
    if WEL_Dict == 0:
        WEL_Dict = {}
    
    # Convert X and Y to Column and Row
    if coordType == 'xy':
        cconvert = lambda x: int(np.floor((x - xll)/cellsize))
        New_WEL[:,0] = np.array([cconvert(xi) for xi in New_WEL[:,0]])
        rconvert = lambda y: int(np.floor((yur - y)/cellsize))
        New_WEL[:,1] = np.array([rconvert(yi) for yi in New_WEL[:,1]])
    
    # If data is in Row/Column format, convert to 0 index
    if coordType == 'rc':
        New_WEL[:,0] = np.array([int(xi) for xi in New_WEL[:,0]-1])
        New_WEL[:,1] = np.array([int(yi) for yi in New_WEL[:,1]-1])

    # Loop through all wells in the dataset to fill dictionary
    for w in range(0,New_WEL.shape[0]):
        r = New_WEL[w,1]
        c = New_WEL[w,0]
                                    
        # Assign flow rate for each well to all stress periods indicated by start and end years
        for per in range(int(New_WEL[w,2]-1),int(New_WEL[w,3]-1)):
            
            try:
                WEL_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult])
            except:
                WEL_Dict[per] = [[LYR,r,c,New_WEL[w,4]*WEL_mult]]
                
    return WEL_Dict

def outputControl(mf):
    # Output Control and Solver
    # Add OC package to the MODFLOW model
    spd = {}
    data2record = ['save head', 'save drawdown', 'save budget', 'print budget']
    for y in range(0,30):
        for m in range(1,13):
            for d in range(0,calendar.monthrange(1984+y,m)[1]):
                spd[y*12 + m - 1, d] = data2record.copy()
#    spd[14,30] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference']
#    for p in [6,20,32,41,54,78,90,102,114,128,138,150,162,175,187,198,213,225,235]:
#        spd[p,0] = data2record.copy()
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
#    oc = flopy.modflow.ModflowOc.load('VM_Test.oc', mf,ext_unit_dict=ext_unit_dict)
    
    # Add PCG package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(mf, mxiter=20, iter1=20)
    
    return oc, pcg


## Build Model
# Assign name and create modflow model object
modelname = 'model\VM_1984'

# Model limits
xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

# Model time
startyear = 1984
endyear = 1985

RCH_PAR = [1.000E-02,5.639E-01,5.000E-01] # Recharge multiplier for urban, natural, and water cover

# Load datasets
ACTIVE_VM_LYR1 = np.loadtxt('ACTIVE_VM_LYR1.asc',skiprows=6) # Active area in layer 1
ACTIVE_VM_LYR2 = np.loadtxt('ACTIVE_VM_LYR2.asc',skiprows=6) # Active area in layer 2
ACTIVE_VM = [ACTIVE_VM_LYR1, ACTIVE_VM_LYR2]

GEO_VM_LYR1 = np.loadtxt('GEO_VM_LYR1.asc',skiprows=6) # Geologic formation distribution in layer 1
GEO_VM_LYR2 = np.loadtxt('GEO_VM_LYR2.asc',skiprows=6) # Geologic formation distribution in layer 2
GEO_VM = [GEO_VM_LYR1, GEO_VM_LYR2]

DEM_VM = np.loadtxt('DEM_VM.asc',skiprows=6) # Digital elevation model

LU = {}
for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
    filename = r'landuse\LU-1984-' + LUtype + '.asc'
    LU[LUtype] = np.loadtxt(filename,skiprows=6)

Precip_Dict = {}
for year in range(int(startyear), int(endyear)):
    for month in range(1,13):
        per = (year - startyear)*12 + month - 1
        filename = r'recharge\Precip_VM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        Precip_Dict[per] = np.loadtxt(filename,skiprows=6)

## Well object
PUMP_array = np.loadtxt(r'PUMP_VM_1984.csv',delimiter=',', skiprows=1, usecols=[0,1,4,5,8]) # pumping in m3 per day


## Organize geologic parameters
formations1 = [0,1] # Number of geologic formations in layer 1
formations2 = [2,3,4,5] # Number of geologic formations in layer 2

ZoneParams = {}
# Zone hydraulic conductivity
ZoneParams['HK'] = [[formations1, formations2],[8.498E+00, 4.320E-04, 2.331E+02, 6.750E-02, 2.518E-01, 8.640E-02]] # m/d

# Zone hydraulic conductivity vertical anisotropy
ZoneParams['VKA'] = [[formations1, formations2],[0.10, 100, 100, 10, 1, 1]]

# Zone specific storage
ZoneParams['SS'] = [[formations1, formations2],[1.000E-06, 6.562E-02, 1.073E-03, 3.176E-02, 1.214E-05, 9.483E-06]] # 1/m

# Zone specific yield
ZoneParams['SY'] = [[formations1, formations2],[1.000E-01, 0.06, 0.15, 0.15, 0.30, 0.01]]

ncol, nrow, mf, dis, bas, lpf = initializeFM(modelname, xll, yll, xur, yur, cellsize, startyear, endyear, ACTIVE_VM, GEO_VM, DEM_VM, ZoneParams)

RCH_DICT = addRecharge(LU_arrays=LU, PRECIP=Precip_Dict, RCH_mult=RCH_PAR, yrzero=startyear, S_YR=startyear, E_YR=endyear)

rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)

WEL_DICT = addNewWells(New_WEL=PUMP_array, LYR=1, coordType='xy', xll=455000, yur=2175000, cellsize=500)

wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT)

oc, pcg = outputControl(mf)

hobs = flopy.modflow.ModflowHob.load('OBS_VM.ob_hob', mf)

# Run Model and post processing
# Write the MODFLOW model input files 
mf.write_input()

# Run the MODFLOW model
success, buff = mf.run_model()
