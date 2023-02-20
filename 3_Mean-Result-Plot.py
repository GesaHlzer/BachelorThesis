# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 18:25:10 2022

_______________________________________________________________________________

These program was made for combining all analysed data of FRAP with THP1 
(beta G-actin citrin) cells and visualise it in a figure
    
It is the 3rd part of 3 single programs:
    - 1_Fitting-Curves-Python.py
    - 2_Analysis_Combine_Resultes.py
    - 3_Mean_Result_Plot.py  
_______________________________________________________________________________

author: Gesa Hölzer
        gesa.hoelzer@t-online.de 
        IAOB, Prof. Eggeling Group
        Friedrich-Schiller-Universität Jena
                  ____________________
                  
This program is for averaging all measured normalized, shifted and photofade corrected 
intensities at each time point,
and plotting the result (intensity means + confidence interval).
Furthermore, an exponential curve ( y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2) ) 
will be plotted as well, using the mean parameters (of y0, A1, τ1, A2, τ2 )
from program 2 (txt-data)
                       --------
(General bleach curve here not plotted, because it was not usable with the acquired data.)
                  ____________________
"""
#-------------------------------------
#Import Packages

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import statistics

SavePlot, MakeBothInOnePlot, SaveBothInOnePlot = False, False, False
#------------------------------------

Loop = ['Control', 'CK666']
tvalue = [2.120, 2.262] # t-Quantil (twosided) for the amount of used curves (here: 9 & 16)

x = 0  # 'Control' data 
x = 1  # 'CK666' data 


Path = 'C:/Users/alleh/OneDrive - Technische Universität Ilmenau/A_UNI/BeachlorArbeit/Python_FRAP_Auswertung/'

# SavePlot = True           # Save plot with either CK-666 or Control sample data

# MakeBothInOnePlot = True    # Plot CK-666 or Control curves in one plot
# SaveBothInOnePlot = True


MeanFitparamFolder = Path + 'UsedData/' + Loop[x] + '/'  # Folder with mean fit parameter file
MicroscopeCurveFolder = Path + 'AdaptedData/' + Loop[x] + '/'  # Folder with normalized, shiftet & photofade corrected data


"""____________________________________________________________________________
Importing adapted microscope curves (CSV) """

# -- read all files with extension ..."PhotofadeCorr.csv" or ..."noCorr.csv" --

usedfiles = glob.glob( MicroscopeCurveFolder + "\*PhotofadeCorr.csv") #Photofade Corrected
# usedfiles = glob.glob( MicroscopeCurveFolder +"\*noCorr.csv")


"""____________________________________________________________________________
Calculate Mean, Median & Standard Derivatian"""


""" (opt.) find minimal acquisition time & define data array size """ 
# ---- use it for importing all data, where the dimension needs to be the same
# --- Part can be skiped if acquisition time was always the same.
# --- Or 'shape' is defined manually

i = 0
tdim = np.zeros(len(usedfiles)) # initialize empty array with its length beeing the number of imported files

for i in range(len(usedfiles)):
    d = (np.loadtxt(usedfiles[i], delimiter = ',', skiprows = 1).shape)
    tdim[i] = d[0] 
    i = i+1

shape = (len(usedfiles),int(min(tdim)),d[1])  # Dimension of AllCurves array


""" Import all 2 dimensional data curves as 3D array """
i = 0
AllCurves = np.zeros(shape)

for i in range(shape[0]):   
    cf_read = np.loadtxt(usedfiles[i], delimiter = ',', skiprows = 1,max_rows=shape[1])
    AllCurves[i,:,:] = cf_read
    i = i+1
    
    
""" mean, std for all """
time = AllCurves[0,:,0] # time pointes
i=0
Mean_recov = np.zeros(shape[1])         # Revovery curve
Std_recov = np.zeros(shape[1])

# Mean_bc = np.zeros(shape[1])          # Bleachcontrol (not useful yet)
# Std_bc = np.zeros(shape[1])

for i in range(shape[1]):
    
    R_recov = statistics.fmean(AllCurves[:,i,1])     # Revovery curve
    Mean_recov[i] = R_recov
    
    S_recov = statistics.stdev(AllCurves[:,i,1])    #	Sample standard deviation
    Std_recov[i] = S_recov
    
    # R_bc = np.mean(AllCurves[:,i,2])  # not implemented
    # Mean_bc[i] = R_two
    # S_bc = np.std(AllCurves[:,i,1])
    # Std_bc[i] = S_one 


"""____________________________________________________________________________
Importing the mean fit parameters (.txt) of the 2nd program and def. the fit"""

meanfitparamfile = open( MeanFitparamFolder + Loop[x] + "-Mean Params of 2-Exp-Fct.txt", "r")

mfp = meanfitparamfile.read().split('\n') # form list of strings to array of floats
mfp = pd.to_numeric(mfp, errors='ignore')


fit = mfp[0] + mfp[1] * (1 - np.exp(- time[3:] /mfp[2])) + mfp[3] * (1 - np.exp(- time[3:] /mfp[4]))


"""____________________________________________________________________________
Create plot of combined results"""

plt.figure(figsize=(10, 5))

plt.scatter(time, Mean_recov, marker='.', label='Mean intensity', color='black')        # Recovery Data

plt.errorbar(time[3:], Mean_recov[3:], yerr = (tvalue[x] * Std_recov[3:]/ np.sqrt(shape[0])), elinewidth = 0.9, linestyle = ' ', color = 'black', capsize=3, capthick=1, label = 'Confidence interval (level: 95%)')
plt.errorbar(time[:3], Mean_recov[:3], yerr = (tvalue[x] * Std_recov[:3]/ np.sqrt(shape[0])), elinewidth = 0.9, color = 'black', capsize=3, capthick=1)

# plt.scatter(time, Mean_bc, marker='.',label='Measurment', color='orange') # Bleachcontrol (not useful yet)

plt.plot(time[3:], fit, color = 'red', label = 'Fit curve with averaged parameters') # fit 


plt.xlabel('Time in s', fontsize=14)
plt.ylabel('Intensity (normalized)', fontsize=16)
plt.grid()
plt.legend(loc = 'best')

if SavePlot == True:
     plt.savefig(Path + 'UsedData/' + Loop[x] + '-Mean-result.png' , format="png")

plt.show()



"""____________________________________________________________________________
_______________________________________________________________________________
Save CK-666 treated and Control in one plot """

if MakeBothInOnePlot == True:
    
    Loop2 = ['CK666' , 'Control'] # CK666 and Control swaped: taking the other

    MeanFitparamFolder2 = Path +'data/' + Loop2[x] + "/"
    MicroscopeCurveFolder2 = Path + 'AdaptedData/' + Loop2[x] + "/" # Folder with 


    """ Importing adapted microscope curves (CSV) """
    """-- read all files with extension ..."PhotofadeCorr.csv" or ..."noCorr.csv" --"""

    usedfiles2 = glob.glob( MicroscopeCurveFolder2 + "\*PhotofadeCorr.csv") # Photofade corrected
    # usedfiles2 = glob.glob( MicroscopeCurveFolder2 +"\*noCorr.csv")


    """ Import data as 3D array """
    i = 0
    AllCurves2 = np.zeros((len(usedfiles2),shape[1],shape[2])) # if error, row amount shorter than in the first part, swap x

    for i in range(len(usedfiles2)):   
        cf_read = np.loadtxt(usedfiles2[i], delimiter = ',', skiprows = 1, max_rows=shape[1])
        AllCurves2[i,:,:] = cf_read
        i = i+1
    
    """ mean, std. for all """
    i=0
    Mean_recov2 = np.zeros(shape[1])   # Revovery curve
    Std_recov2 = np.zeros(shape[1])

    for i in range(shape[1]):
    
        R_recov2 = statistics.fmean(AllCurves2[:,i,1])
        Mean_recov2[i] = R_recov2
    
        S_recov2 = statistics.stdev(AllCurves2[:,i,1])
        Std_recov2[i] = S_recov2
    

    """ Importing the mean fit parameters (.txt) of the 2nd program and def. the fit """

    meanfitparamfile = open( MeanFitparamFolder2 + "Mean Params of 2-Exp-Fct.txt", "r")

    mfp2 = meanfitparamfile.read().split('\n')
    mfp2 = pd.to_numeric(mfp2, errors='ignore')
    mfp2 = np.delete(mfp2, 5)

    fit2 = mfp2[0] + mfp2[1] * (1 - np.exp(- time[3:] /mfp2[2])) + mfp2[3] * (1 - np.exp(- time[3:] /mfp2[4]))


    """________________________________________________________________________
    Create plot of combined results"""

    plt.figure(figsize=(10, 5))
    
    plt.scatter(time, Mean_recov, marker='.', label='Mean intensities of '+Loop[x], color='black')    # Recovery Data   
    plt.scatter(time, Mean_recov2, marker='.', label='Mean intensities of '+Loop2[x], color='navy')   
    
    plt.errorbar(time[3:], Mean_recov[3:], yerr = tvalue[x] * Std_recov[3:]/ np.sqrt(shape[0]), elinewidth = 0.5, linestyle = ' ', color = 'black', capsize=3, capthick=1, label = '95 % confidence interval (' + Loop2[x] + ')')
    plt.errorbar(time[:3], Mean_recov[:3], yerr = tvalue[x] * Std_recov[:3]/ np.sqrt(shape[0]), elinewidth = 0.5, color = 'black', capsize=3, capthick=1)
   
    plt.errorbar(time[3:], Mean_recov2[3:], yerr = tvalue[x-1] * Std_recov2[3:]/ np.sqrt(len(usedfiles2)), elinewidth = 0.5, linestyle = ' ', color='navy', capsize=3, capthick=1, label = '95% confidence interval (' + Loop2[x] + ')')
    plt.errorbar(time[:3], Mean_recov2[:3], yerr = tvalue[x-1] * Std_recov2[:3]/ np.sqrt(len(usedfiles2)), elinewidth = 0.5, color='navy', capsize=3, capthick=1)
    
    plt.plot(time[3:], fit, color = 'red', label= 'Fit with averaged parameters (' +Loop[x]+')') # Fit 
    plt.plot(time[3:], fit2, color = 'crimson', label = 'Fit with averaged parameters (' +Loop2[x]+')') # Fit 

    plt.legend(loc = 'best')
    plt.xlabel('Time in s', fontsize=14)
    plt.ylabel('Intensity (normalized)', fontsize=16)
    plt.grid()

    if SaveBothInOnePlot == True:
        plt.savefig(Path + 'UsedData/' + 'Mean-result.png' , format="png")
    plt.show()


###############################################################################

print('___________')
print('done')




