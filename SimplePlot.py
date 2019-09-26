# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:02:08 2018

@author: MM
"""
import flopy
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import calendar
import seaborn as sns
from pathlib import Path
import os

# Set this to the path where the file is located
filelocation = Path(r'D:\MMautner\flopyTutorial')

s_name = 'VM_1984'

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

STRT_YEAR = 1984
END_YEAR = 2014

ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

# Load datasets
ACTIVE_VM_LYR1 = np.loadtxt(filelocation / 'ACTIVE_VM_LYR1.asc',skiprows=6)
ACTIVE_VM_LYR2 = np.loadtxt(filelocation / 'ACTIVE_VM_LYR2.asc',skiprows=6)
DEM_VM = np.loadtxt(filelocation / 'DEM_VM.asc',skiprows=6)
GEO_VM = np.loadtxt(filelocation / 'GEO_VM_LYR2.asc',skiprows=6)

# Head Dictionary
hds = bf.HeadFile(filelocation / 'model' / (s_name + '.hds'))
    
#%% Heads Contour
fig, axes = plt.subplots(2, 2, figsize=(7,6.3))
plt.set_cmap('rainbow_r')
axes = axes.flat
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])

mapTitle = ['March','June','September','December']

# Choose time-step and period for 4 different contour maps: March 31, June 30, September 30, and December 31
mapkstpkper = [(30,2),(29,5),(29,8),(30,11)]

for a, subtitle in enumerate(mapTitle):
    h = hds.get_data(mflay=1, kstpkper=mapkstpkper[a])
    h[ACTIVE_VM_LYR2!=1] = np.nan
    
    im = axes[a].imshow(h)
    CS = axes[a].contour(ACTIVE_VM_LYR1, colors='k', linewidths=2)
    axes[a].xaxis.set_visible(False)
    axes[a].yaxis.set_visible(False)
    axes[a].set_title(mapTitle[a].format(a+1))
    
fig.subplots_adjust(right=0.8)
fig.colorbar(im, cax=cbar_ax, label='Groundwater Head (m)')

#%% Budget
mf_list = flopy.utils.MfListBudget(filelocation / 'model' / (s_name + '.list'))
incremental, cumulative = mf_list.get_budget()

df_1Bdget, df_extra = mf_list.get_dataframes(start_datetime='12-31-1983')
    
mthly_Bdget = df_1Bdget.drop(['CONSTANT_HEAD_IN','TOTAL_IN','CONSTANT_HEAD_OUT','RECHARGE_OUT','TOTAL_OUT','IN-OUT','PERCENT_DISCREPANCY','WELLS_IN'], axis=1)

mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
mthly_Bdget = mthly_Bdget.multiply(30/1000000)

cols = mthly_Bdget.columns.tolist()
# reorder columns
cols = [cols[1]] + [cols[0]] + [cols[3]] + [cols[2]] 
# "commit" the reordering
mthly_Bdget = mthly_Bdget[cols]

ax = mthly_Bdget['04-01-1984':'12-31-1985'].plot.area(stacked=True,figsize=(8,9),color=['blue','lightblue','red','lightpink'])

plt.ylabel(r'Volume ($hm^3$)')
plt.title('Groundwater Budget')
plt.legend(['Precipitation','Storage: In','Pumping','Storage: Out'],loc=4)

#%% Cumulative overdraft
plt.figure(figsize=(12,7.2))

i = 0

# Scale datasets by 1 million m3
df_extra['IN'] = df_extra['RECHARGE_IN'].divide(1000000)
df_extra['OUT'] = df_extra['WELLS_OUT'].divide(1000000)

# Take the difference between recharge and pumping
df_extra['INOUTCUMSUM'] = df_extra['IN'] - df_extra['OUT']

df_extra.INOUTCUMSUM['01-01-1984':'12-31-1984'].plot()

plt.ylabel(r'Volume (million m$^3$)')
plt.title('Cumulative In - Out')

#%% Time series by location
t = pd.DatetimeIndex(freq='D',start='01/01/1984',end='12/31/1984')

# Choose row and column of hydrograph location
coords = [29,77] # pump, [52,61] # subs, [90,25] # mtn, [59,54] #

# Create an empty vector
hTimeS = np.zeros(366)

# Loop through time steps in each stress period
j = 0
for y in range(0,1):
    for m in range(0,12):
        for d in range(0,(calendar.monthrange(1984+y,m+1)[1])):
            h = hds.get_data(kstpkper=(d,y*12+m),mflay=1)
            hTimeS[j] = h[coords[0],coords[1]]
            j+=1

plt.plot(t, hTimeS)
    
plt.xlabel('Date')
plt.ylabel('Head Elevation')

#%% Head Surface

# Initialize figure
fig = plt.figure(figsize=(12,7))
ax = fig.gca(projection='3d')

# Extract heads data from Layer 1 and Layer 2 (zero indexed as mflay = 0 and 1). Data is for last time step: stress period 11, time-step 31
h = hds.get_data(mflay=[0,1],kstpkper=(30,11))

# Create a matrix with all values equal to the smallest value above 0, this is to filter out dry cells that will be a very negative number due to MODFLOW conventions
HNEW = np.ones(h.shape)*min(h[h>0])

# Fill in values for each layer that are not dry cells
for i, head in enumerate(h):
    HNEW[i][head>-900] = head[head>-900]

# plot Layer 2 as a 3D wireframe like in the example mplot3d/wire3d_demo, change to HNEW[0] for Layer 1
Z = HNEW[1]
x = np.arange(0,ncol*cellsize,cellsize)
y = np.arange(0,nrow*cellsize,cellsize)
X, Y = np.meshgrid(x, y)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.RdYlBu, vmin=2000, vmax=3000)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5, label='Elevation (masl)')
plt.xlabel('Easting')
plt.ylabel('Northing')