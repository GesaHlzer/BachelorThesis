# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 09:2:15 2022

_______________________________________________________________________________

These program was made for the data analysis of single measurements 
    of FRAP with THP1 (beta G-actin citrin) cells
    
It is the 2nd part of 3 single programs:
    - 1_Fitting-Curves-Python.py
    - 2_Analysis_Combine_Resultes.py
    - 3_Mean_Result_Plot.py 
_______________________________________________________________________________

author: Gesa Hölzer
        gesa.hoelzer@t-online.de 
        IAOB, Prof. Eggeling Group
        Friedrich-Schiller-Universität Jena
                  ____________________
    
These program is for the analysis of the recovery-curves,
in order to do that you import csv-data, fit the curve, decide the part of data you want to include and what you want to save.
It automatically makes a folder with your results in it.
For using program 2 you need to save the analysed data.
        
                  
Define in this program :

     - your measurement folder (file_path) & the file you want to analyse (file_name)
     - If you want to exlude data points of the end of the measurement (m)
     - What you want to plot/ save
     - If you take a bleach-control into account (and from to where to fit the bleach function: x1, x2)
                  ____________________
        
FRAP measurement standards in order to use this program as it is:
    
    - Image min. 60 s recovery (~200 cycles)
    - Acquisition interval 0.3 s 
    - Time unit in s (not ms)!
    - 3 cycles before bleaching
    - R1 Bbleach acquisition, R2 bleachcontrol acquisition
                  ____________________

"""

# -------------------------------------
# Import Packages
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as spo
from scipy.optimize import curve_fit
import os.path
import glob
# -------------------------------------

#Inititialise variables for later
TwoExp, SaveData, SaveTwoExp, CombinedPlot, SaveCombindedPlot, OneExp, SaveOneExp, ExpA, SaveExpA, Bleachcontrol, IncludeBleachcontrol, CopyMicroscopeData  = False, False, False, False, False, False, False, False, False, False, False, False
x, f = 0, 1

"""___________________________________________________________________________________________________________________________
__________________________________THIS BLOCK NEEDS TO BE CUSTOMISED!__________________________________________________________
______________________________________________________________________________________________________________________________"""

TimeAxisIn_ms = False # If the aquisition intervall was mesured with a 300 ms intervall instead of 0.3 s, set to True

Path = 'C:/Users/alleh/OneDrive - Technische Universität Ilmenau/A_UNI/BeachlorArbeit/Python_FRAP_Auswertung/'
Folder = 'UsedData/'

# --------- Folder with control or CK-666 treated measurement data --------
Loop = ['Control', 'CK666']
x = 0   # 'Control' data folder
x = 1   # 'CK666' data folder

f = 16  # Choose number of file in folder (f_min = 1) 

m = 1  # amount of data points (m > 0) excluded from fitting in the end of the curve 

w = 2   # customise x_axis length in Plots, it will be the seconds you measured + w sec. (standard: 2)

""" --------- Uncomment which fits and plots you want to show and save (Ctr+1 or remove "# ") --------------"""

TwoExp = True       #required to start adjusting the fit

SaveTwoExp = True

SaveData = True           # Save data (fitted & calculated) in csv
CopyMicroscopeData = True # Copy used microscope csv file to an extra folder 
#                             #(for programm 3, remember which measurement files have been used)

# OneExp = True             # Show plot of exp. fit 1st Order
# SaveOneExp = True 

# ExpA = True               # Show plot of exp. fit to power of alpha
# SaveExpA = True

CombinedPlot = True       # Show plot of all three fits in one figure
SaveCombindedPlot = True  

""" Bleachcontrol """
Bleachcontrol = True        # if not mesured (as R2): out comment this line (#)
IncludeBleachcontrol = True # if it shall not be included: out comment this line (#)

x1 = 78     # amount of data points to exclude in the beginning; x1 > 0
x2 = 5      # amount of data points to exclude in the end; standard: m


"""________________________________________________________________________________________________________________________
____________________________DONE WITH CUSTOMISATION________________________________________________________________________
___________________________________________________________________________________________________________________________"""

# --- Initialising saving folders ---

ResultFolder = Path + Folder + Loop[x]  + '/Result Data/'  # Where to save the result parameters
PlotFolder = Path + Folder + Loop[x] + '/Result Plots/'  # Where to save the plots

""" ___________________________________________________________________________
Include the by 'f' determined file 
"""

# --- Creating a list with the names of all csv data in the folder ---
filenames = glob.glob(Path + Folder + Loop[x] + "/" + "\*.csv")


f = f-1     # Start counting from 1 (instead 0)

file = filenames[f]
file_name = file[(len(Path)+len(Folder)+len(Loop[x])+1):]

z = 0       # Initialize counter
for z in range(len(filenames)):
    print_file = filenames[z]
    print_file = print_file[(len(Path)+len(Folder)+len(Loop[x])+1):]
    print(print_file)
    
# Print information for user  
print('')
print('Amount of files in folder:', len(filenames))
print('')
print('use for f an int between 1 and ', len(filenames))
print('')
print('Using file: ', file_name)
print('')


"""____________________________________________________________________________
Import relevant part of the LSM980 csv-table as an array
"""

if Bleachcontrol == True or IncludeBleachcontrol == True:
    df = pd.read_csv(file, delimiter=",", usecols=[0, 5, 8], skiprows=1)
else:
    df = pd.read_csv(file, delimiter=",", usecols=[0, 5], skiprows=1)

data = df.to_numpy() # Convert the entire DataFrame tu numpy array

if TimeAxisIn_ms == True:  # If the aquisition intervall was mesured with a 300 ms intervall instead of 0.3 s
    data[:, 0] = data[:, 0]/1000
    

data_length = len(data[:, 0]) # Amount of data points (with 3 bleaching rows)
# print(data_length)


"""____________________________________________________________________________
Implement of bleach funktion (from an extra measurement) 
"""

if Bleachcontrol == True:

    def linear_bleach(t, y, b):
        return y + b * t            # Define fit function 
    
    popt_lin_bl, pcov = spo.curve_fit(linear_bleach, data[x1:-x2, 0], data[x1:-x2, 2], p0=[data[0, 2], 1], method='lm')   # Fit, lm=levenberg-marquardt
    popt_lin_bl
    # print(popt_lin_bl)

    bleach_relation = linear_bleach(data[:, 0], *popt_lin_bl) / data[x1, 2] # Calculate photofade factor for each time point (photofading factor)
    bleach_correction = 1/bleach_relation 

    ndata = np.vstack((data[:, 0], np.multiply(data[:, 1], bleach_correction[:]))).transpose() # Adapt data
    
    
    
    """ --- Plot recovery data, bleachcontrol data and bleachimplement in a graph ---"""
    
    plt.figure(figsize=(10, 5))
    
    plt.scatter(data[:, 0], data[:, 1], marker='.', # Recovery Data 
                label='Measurment', color='blue') 
    plt.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                label='Measurment', color='black') # Bleachimplement
    plt.scatter(data[:, 0], data[:, 2], marker='.',
                label='Measurment', color='orange') # Bleachcontrol
    plt.plot(data[x1:-x2, 0], linear_bleach(data[x1:-x2, 0], *popt_lin_bl)) # Linear fit of bleachcontrol
    
    plt.xlabel('Time in s')
    plt.grid()
    
    if not os.path.isdir(PlotFolder): # If required, create result plot folder
        os.makedirs(PlotFolder)
        
    plt.savefig(PlotFolder+file_name[:-4]+'-Bleachimplement.png', format="png") # Save plot
    plt.show()
   

    if IncludeBleachcontrol == False: # No adaptation
        ndata = data[:, :2]
else:
    ndata = data # No adaptation

"""____________________________________________________________________________
Translate the start of the recovery curve (IA) to the (0,0)-coordinate 
"""
ndata = ndata[:, :] - ndata[3, :]

"""____________________________________________________________________________
Normalisation (Min-Max)
"""
norm = (ndata[0, 1] + ndata[1, 1] + ndata[2, 1])/3 # Averaging the first 3 measured values

ndata = np.hstack((ndata[:, 0], ndata[:, 1]/norm))
ndata = ndata.reshape(2, data_length).transpose() 


"""____________________________________________________________________________
 Save (Copy) adapted and used data to a single, additional folder (overfolder) 
 """
 
if CopyMicroscopeData == True:
    
    dn = pd.DataFrame(data=ndata)
    
    if not os.path.isdir(Path + '/AdaptedData/' + Loop[x] + "/" ):
        os.makedirs(Path + '/AdaptedData/'+ Loop[x] + "/" )
    
    if IncludeBleachcontrol == True:
        dn.to_csv(Path + '/AdaptedData/' + Loop[x] + "/" + file_name[:-4]+'+PhotofadeCorr.csv', index=None)
    else:
        dn.to_csv(Path + '/AdaptedData/' + Loop[x] + "/" + file_name[:-4]+'+noCorr.csv', index=None)


"""____________________________________________________________________________
Fit of the 2-Component Exponential Function 
"""
m = -m # "-" makes the counter start backwards from the last measured value

# Recovery part of data
t = ndata[3:m, 0]  # Time
F = ndata[3:m, 1]  # R1: bleached part
# print(F)

# --------- Define Fit Function -------------------

def exp_decay_two(t, y_2, A_1, t_1, A_2, t_2):
    return y_2 + A_1 * (1 - np.exp(-t / t_1)) + A_2 * (1 - np.exp(-t / t_2))


# ------------ Set Starting Parameters ------------------------

guess_two = [0, 23, 4, 9, 26]
#         = [y_2, A_1,t_1,A_2,t_2]

# ------------------- Fitting ---------------------------

popt_two, pcov = spo.curve_fit(exp_decay_two, t, F, p0=guess_two, method='lm')
popt_two
d_pop_two = np.sqrt(np.diag(pcov))

y_2, A_1, t_1, A_2, t_2 = popt_two
d_y_2, d_A_1, d_t_1, d_A_2, d_t_2 = d_pop_two

if TwoExp == True:

    print()
    print('For the fit function y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)):')
    print()
    print('y0=', y_2, ' +- ', d_y_2)
    print('A1=', A_1, ' +- ', d_A_1)
    print('τ1=', t_1, ' +- ', d_t_1)
    print('A2=', A_2, ' +- ', d_A_2)
    print('τ2=', t_2, ' +- ', d_t_2)
    print()

# ------------ Calculate Residuals -----------------

Residual = abs(F - exp_decay_two(t, *popt_two)) # Difference between each measurement value and fit (in absolute values)
Res_Mean = np.mean(Residual)

"""____________________________________________________________________________
Calculate adiddional parameters:
"""
IB = 1
IA = ndata[3, 1]  # = 0 = y_0
# IB = 1 included with normalisation

IE = A_1 + A_2

k_1 = 1/t_1
k_2 = 1/t_2
t_1_half = - math.log(0.5)*t_1
t_2_half = - math.log(0.5)*t_2

# I_mobile = (IE - IA)/(IB - IA) = IE
I_immobile = 1 - IE

# Fractions 
F_1 = A_1/IE 
F_2 = A_2/IE 
F_immobile = I_immobile 
F_mobile = IE 

if TwoExp == True:
    print('IE=', IE)
    print('I_immobile=', I_immobile)
    print('k1=', k_1)
    print('τ1_half=', t_1_half)
    print('k2=', k_2)
    print('τ2_half=', t_2_half)
    print('F1=', F_1)
    print('F2=', F_2)
    print()


"""____________________________________________________________________________
Plot data & fit funktion
"""
# 2
if TwoExp == True:
    
    # Plot (ax) and Residuals (Res)
    fig2, (ax, Res) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, figsize=(10, 7), constrained_layout=True) 
    fig2.suptitle('Fit with Two Exponential Summands', fontsize=20)

    if IncludeBleachcontrol == True:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement (photofade corrected)', color='k')
    else:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement', color='k')

    ax.plot(t, exp_decay_two(t, *popt_two), 'g-',
            label='Fit: y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)) \n y0=%5.3f, A1=%5.3f, τ1=%5.3f,\n A2=%5.3f, τ2=%5.3f' % tuple(popt_two))
    ax.grid()
    ax.set_ylabel('Intensity (normalized)', fontsize=16)
    #ax.set_xlabel('Time in s', fontsize=14)
    ax.set_xlim(-5, ndata[-1, 0]+w)
    ax.legend(loc='best', fontsize=14)
    

    Res.plot(t, Residual)
    Res.set_ylabel('Residuals', fontsize=16)
    Res.set_xlabel('Time in s', fontsize=16)
    Res.axhline(y=Res_Mean, xmin=0, xmax=60, color='r', linewidth=0.5,
                label='Mean difference: %5.3f' % float(Res_Mean))
    Res.legend(loc='best', fontsize=14)
    Res.set_xlim(-5, ndata[-1, 0]+w)
    Res.grid()

    fig2.align_ylabels((ax, Res)[:])
    fig2.tight_layout()
    fig2.show()

    if SaveTwoExp == True:
        
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
            
        if IncludeBleachcontrol ==  True:
            fig2.savefig(PlotFolder+file_name[:-4]+'---Fit-2nd-Order-Plot_+PhotofadeCorr.pdf', format="pdf")
            
        else: # different file name
            fig2.savefig(PlotFolder+file_name[:-4]+'---Fit-2nd-Order-Plot_noCorr.pdf', format="pdf")


"""____________________________________________________________________________
Opt.: Evaluate Other Fit Funktions 
"""

# ------------------- Fit with 1 exponential function-----------------------------------

# Define Fit Function

def exp_decay_one(t, y_1, A_0, t_0):
    return y_1 + A_0 * (1 - np.exp(-t / t_0))

# Set Starting Parameters

guess_one = [0, 22, 5]
#         = [ B_1, A_0, t_0]

# Fit

popt_one, pcov = curve_fit(exp_decay_one, t, F, p0=guess_one, method='lm') 
popt_one
d_popt_one = np.sqrt(np.diag(pcov)) # Standard deviation

y_1, A_0, t_0 = popt_one
d_y_1, d_A_0, d_t_0 = d_popt_one

if OneExp == True:
    print('For the fit function y0 +  A * (1 - exp(-t/τ)) :')
    print('y0=', y_1, ' +- ', d_y_1)
    print('A=', A_0, ' +- ', d_A_0)
    print('τ=', t_0, ' +- ', d_t_0)
    print()

# Residuals
Residual_one = abs(F - exp_decay_one(t, *popt_one))
Res_Mean_one = np.mean(Residual_one)

# Create Plot
if OneExp == True:

    fig1, (ax, Res) = plt.subplots(2, 1, gridspec_kw={
        'height_ratios': [4, 1]}, figsize=(10, 7), constrained_layout=True)
    fig1.suptitle('Fit with One Exponential Summand', fontsize=20)
    
    if IncludeBleachcontrol == True:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.', label='Measurement (photofade corrected)', color='k')
    else:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.', label='Measurement', color='k')
        
    ax.plot(t, exp_decay_one(t, *popt_one), 'b-', label='Fit: y0 + A * (1 - exp(-t/τ))\n y0=%5.3f, A=%5.3f,τ=%5.3f' % tuple(popt_one))
    ax.grid()
    ax.set_ylabel('Intensity (normalized)', fontsize=16)
    # ax.set_xlabel('Time in s', fontsize=14)
    ax.set_xlim(-5, ndata[-1, 0]+w)
    ax.legend(loc='best', fontsize=14)

    Res.plot(t, Residual_one)
    Res.set_ylabel('Residuals', fontsize=16)
    Res.set_xlabel('Time in s', fontsize=16)
    Res.axhline(y=Res_Mean_one, xmin=0, xmax=60, color='r', linewidth=0.5,
                label='Mean difference: %5.3f' % float(Res_Mean_one))
    Res.legend(loc='best', fontsize=14)
    Res.set_xlim(-5, ndata[-1, 0]+w)
    Res.grid()

    fig1.align_ylabels((ax, Res)[:])
    fig1.tight_layout()
    fig1.show()

    if SaveOneExp == True:
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
            
        if IncludeBleachcontrol ==  True:
            fig1.savefig(PlotFolder+file_name[:-4]+'---Fit-1st-Order-Plot_+PhotofadeCorr.pdf', format="pdf")
        else:
            fig1.savefig(PlotFolder+file_name[:-4]+'---Fit-1st-Order-Plot_noCorr.pdf', format="pdf")


# -------------------------- One exponatial summand to the power of alpha --------------------------------------

# Define Fit Function

def exp_decay_a(t, y_a, A_a, t_a, alpha):
    return y_a + A_a * (1 - np.exp(-(t/t_a)**alpha)) # without setting it here: alpha < 1

# Set Starting Parameters
guess_a = [0, 30, 6, 0.7]
#       = [y_a, A_a, t_a, alpha]

# Fit
popt_a, pcov = curve_fit(exp_decay_a, t, F, p0=guess_a, bounds=(-5, 500))
popt_a
d_popt_a = np.sqrt(np.diag(pcov))
y_a, A_a, t_a, alpha = popt_a
d_y_a, d_A_a, d_t_a, d_alpha = d_popt_a

if ExpA == True:
    print('For the fit function y0 + A * (1 - exp((-t/τ)**a) ) :')
    print('y0=', y_a, ' +- ', d_y_a)
    print('A=', A_a, ' +- ', d_A_a)
    print('τ=', t_a, ' +- ', d_t_a)
    print('alpha =', alpha, ' +- ', d_alpha)
    print()

# Residuals
Residual_a = abs(F - exp_decay_a(t, *popt_a))
Res_Mean_a = np.mean(Residual_a)

if ExpA == True:
    figA, (ax, Res) = plt.subplots(2, 1, gridspec_kw={
        'height_ratios': [4, 1]}, figsize=(10, 7), constrained_layout=True)
    figA.suptitle('Fit with One Exponential Summand to Power of α', fontsize=20)

    if IncludeBleachcontrol == True:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement (photofade corrected)', color='k')
    else:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement', color='k')

    ax.plot(t, exp_decay_a(t, *popt_a), 'y-',
            label='Fit of y0 + A * (1 - exp((-t/τ)**α)): \n y0=%5.3f, A=%5.3f, τ=%5.3f, α=%5.3f' % tuple(popt_a))
    ax.grid()
    ax.set_ylabel('Intensity (normalized)', fontsize=16)
    # ax.set_xlabel('Time in s', fontsize=14)
    ax.set_xlim(-5, ndata[-1, 0]+w)
    ax.legend(loc='best', fontsize=14)

    Res.plot(t, Residual_a)
    Res.set_ylabel('Residuals', fontsize=16)
    Res.axhline(y=Res_Mean_a, xmin=0, xmax=60, color='r', linewidth=0.5,
                label='Mean difference: %5.3f' % float(Res_Mean_a))
    Res.legend(loc='best', fontsize=14)
    Res.set_xlabel('Time in s', fontsize=16)
    Res.set_xlim(-5, ndata[-1, 0]+w)
    Res.grid()

    figA.align_ylabels((ax, Res)[:])
    figA.tight_layout()
    figA.show()

    if SaveExpA == True:
        
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
            
        if IncludeBleachcontrol ==  True:
            figA.savefig(PlotFolder+file_name[:-4]+'---Fit-Exp-Alpha-Plot_+PhotofadeCorr.pdf', format="pdf")
        else:
            figA.savefig(PlotFolder+file_name[:-4]+'---Fit-Exp-Alpha-Plot_noCorr.pdf', format="pdf")

"""____________________________________________________________________________
Combined Plot
"""

if CombinedPlot == True:

    figC, ([ax1, ax2, ax3], [Res1, Res2, Res3]) = plt.subplots(2, 3, figsize=(
        30, 7), gridspec_kw={'height_ratios': [4, 1]}, layout='constrained')
    figC.suptitle('Comparison of Different Fitting Functions', fontsize=20)

    # ______axes_______
    ax1.set_ylabel('Intensity (normalized)', fontsize=16)
    ax2.set_ylabel('Intensity (normalized)', fontsize=16)
    ax3.set_ylabel('Intensity (normalized)', fontsize=16)
    Res1.set_ylabel('Residuals', fontsize=16)
    Res2.set_ylabel('Residuals', fontsize=16)
    Res3.set_ylabel('Residuals', fontsize=16)
    Res1.set_xlabel('Time in s', fontsize=16)
    Res2.set_xlabel('Time in s', fontsize=16)
    Res3.set_xlabel('Time in s', fontsize=16)
    ax1.set_xlim(-5, ndata[-1, 0]+w)
    ax2.set_xlim(-5, ndata[-1, 0]+w)
    ax3.set_xlim(-5, ndata[-1, 0]+w)
    Res1.set_xlim(-5, ndata[-1, 0]+w)
    Res2.set_xlim(-5, ndata[-1, 0]+w)
    Res3.set_xlim(-5, ndata[-1, 0]+w)

    figC.align_ylabels((ax1, Res1)[:])

    # ______1______
    ax1.set_title('One Exponential Summand')
    ax1.plot(t, exp_decay_one(t, *popt_one), 'b-',
             label='Fit: y0 + A * (1 - exp(-t/τ)): \n y0=%5.3f, A=%5.3f,τ=%5.3f' % tuple(popt_one))
    if IncludeBleachcontrol == True:
        ax1.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement (photofade corrected)', color='k')
    else:
        ax1.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement', color='k')
    ax1.grid(True)
    ax1.legend(loc='best', fontsize=14)

    Res1.plot(t, Residual_one)
    Res1.set_xlabel('Time in s')
    Res1.grid(True)
    Res1.axhline(y=Res_Mean_one, xmin=0, xmax=60, color='r', linewidth=0.5,
                 label='Mean difference: %5.4f' % float(Res_Mean_one))
    Res1.legend(loc='best', fontsize=14)

    # ____2_______
    ax2.set_title('Two Exponential Summands')
    ax2.plot(t, exp_decay_two(t, *popt_two), 'g-',
             label='Fit: y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)) : \n y0=%5.3f, A1=%5.3f, \n τ1=%5.3f, A2=%5.3f, τ2=%5.3f' % tuple(popt_two))  # fit
    if IncludeBleachcontrol == True:
        ax2.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement (photofade corrected)', color='k')
    else:
        ax2.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement', color='k')
    ax2.grid(True)
    ax2.legend(loc='best', fontsize=14)

    Res2.plot(t, Residual)
    Res2.set_xlabel('Time in s')
    Res2.axhline(y=Res_Mean, xmin=0, xmax=60, color='r', linewidth=0.5,
                 label='Mean difference: %5.4f' % float(Res_Mean))
    Res2.grid(True)
    Res2.legend(loc='best', fontsize=14)

    # _____alpha______
    ax3.set_title('Exponential Function to Power of α')
    ax3.plot(t, exp_decay_a(t, *popt_a), 'y-',
             label='Fit: y0 + A * (1 - exp((-t/τ)**α)): \n y0=%5.3f, A=%5.3f, \n τ=%5.3f, α=%5.3f' % tuple(popt_a))
    if IncludeBleachcontrol == True:
        ax3.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement (photofade corrected)', color='k')
    else:
        ax3.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement', color='k')
    ax3.grid(True)
    ax3.legend(loc='best', fontsize=14)

    Res3.plot(t, Residual_a)
    Res3.set_xlabel('Time in s')
    Res3.grid(True)
    Res3.axhline(y=Res_Mean_a, xmin=0, xmax=60, color='r', linewidth=0.5,
                 label='Mean difference: %5.4f' % float(Res_Mean_a))
    Res3.legend(loc='best', fontsize=14)

    figC.show()

    if SaveCombindedPlot == True:
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
        if IncludeBleachcontrol ==  True:
            figC.savefig(
                PlotFolder+file_name[:-4]+'---Fit-Plots-Combined_+PhotofadeCorr.pdf', format="pdf")
        else:
            figC.savefig(
                PlotFolder+file_name[:-4]+'---Fit-Plots-Combined_noCorr.pdf', format="pdf")

"""____________________________________________________________________________
Save Parameters 
"""
if SaveData == True:

    # # As List
    # Result = [['Exp.function 2nd order', 0], #0
    #           ['IE', IE],
    #           ['I_immobile', I_immobile],
    #           ['y0', y_2],
    #           ['F1 [%]', F_1],
    #           ['A1', A_1],
    #           ['t1', t_1],
    #           ['t1_half', t_1_half],
    #           ['k1', k_1],
    #           ['F2 [%]', F_2],
    #           ['A2', A_2],
    #           ['t2', t_2],
    #           ['t2_half', t_2_half],
    #           ['k2', k_2],
    #           ['F_mobile [%]', F_mobile],
    #           ['F_immobile [%]', F_immobile], #15
    #           ['----------',0],
    #           ['Opt.: y0 +  A * (1 - exp(-t/t0))', 0],
    #           ['y0', y_1],
    #           ['A', A_0],
    #           ['t0', t_0],
    #           ['---------',0],
    #           ['Opt.: y0 + A * (1 - exp((-t/t0)**a))', 0],
    #           ['y0', y_a],
    #           ['A', A_a],
    #           ['t0', t_a],
    #           ['alpha', alpha]]

    # As  Numpy Array
    Result = np.array([('Exp. function 2nd order', 0),  # 0
                       ('IE', IE),
                       ('I_immobile', I_immobile),
                       ('y0', y_2),
                       ('F1 [%]', F_1),
                       ('A1', A_1),
                       ('t1', t_1),
                       ('t1_half', t_1_half),
                       ('k1', k_1),
                       ('F2 [%]', F_2),
                       ('A2', A_2),
                       ('t2', t_2),
                       ('t2_half', t_2_half),
                       ('k2', k_2),
                       ('F_mobile [%]', F_mobile),
                       ('F_immobile [%]', F_immobile),  # 15
                       (' ', 0),
                       ('Opt.: y0 +  A * (1 - exp(-t/t0))', 0),
                       ('y0', y_1),
                       ('A', A_0),
                       ('t0', t_0),
                       (' ', 0),
                       ('Opt.: y0 + A * (1 - exp((-t/t0)**a))', 0),
                       ('y0', y_a),
                       ('A', A_a),
                       ('t1', t_a),
                       ('alpha', alpha),
                       (' ', 0),
                       ('m', -m),
                       ('x1', x1),
                       ('x2', x2)],
                      dtype=('a8, f4'))

    # ResultFolder    
    if not os.path.isdir(ResultFolder):
        os.makedirs(ResultFolder)    
    if IncludeBleachcontrol == True:
        dR = pd.DataFrame(data=Result)
        dR.to_csv(ResultFolder+'Result-Parameters_' +
                  file_name[:-4]+'+PhotofadeCorr.csv', index=None)
    else:
        dR = pd.DataFrame(data=Result)
        dR.to_csv(ResultFolder+'Result-Parameters_'+file_name, index=None)


"""____________________________________________________________________________
Dictonary with all Results
"""

ResultDict = {'Exp.function 2nd order': {'IE': IE,
                                         'I_immobile': I_immobile,
                                         'y0': [y_2, d_y_2],
                                         'F1 [%]': F_1,
                                         'A1': [A_1, d_A_1],
                                         't1': [t_1, d_t_1],
                                         't1_half': t_1_half,
                                         'k1': k_1,
                                         'F2 [%]': F_2,
                                         'A2': [A_2, d_A_2],
                                         't2': [t_2, d_t_2],
                                         't2_half': t_2_half,
                                         'k2': k_2,
                                         'F_mobile [%]': F_mobile,
                                         'F_immobile [%]': F_immobile},
              'Exp.function 1st order': {'y0': [y_1, d_y_1],
                                         'A1': [A_0, d_A_0],
                                         't1': [t_0, d_t_0]},
              'Exp.function with power to a': {'y0': [y_a, d_y_a],
                                               'A': [A_a, d_A_a],
                                               't0': [t_a, d_t_a],
                                               'α': [alpha, d_alpha]}
              }

##############################################################################

print('___________')
print('done')
