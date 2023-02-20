# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:06:20 2022

_______________________________________________________________________________

This program was made for the analysation of FRAP curves with 
    THP1 (beta G-actin citrin) cells
    
It is the 2nd of 3 programs:
    - 1_Fitting-Curves-Python.py
    - 2_Analysis_Combine_Resultes.py
    - 3_Mean_Result_Plot.py 
_______________________________________________________________________________

author: Gesa Hölzer
        gesa.hoelzer@t-online.de OR gesa.hoelzer@uni-jena.de
        IAOB, Prof. Eggeling Group
        Friedrich-Schiller-Universität Jena
                  ____________________
       
If not already in this format: Save all control and CK-666 treated data direct in different folders and name them
'Control' & 'CK666'

   # or change folder names in line 34 in Loop (was initialised to control which data you analyze -> change x)    
                 _____________________
   
This program is for calculating & saving automaticly 
- the mean, standard derivation & median of the imported data (from Program 1) as an Excel sheet
- the mean of the pararameters of the two-exponatial-term Curve as txt (for Program 3)

                  ____________________
                  
"""
import numpy as np
import pandas as pd
import glob
import statistics

"""________________Customize_x_______________________________________________
 ____________________________________(Which data you want?) ___________________"""
# getting csv files from the folder MyProject

Path = 'C:/Users/alleh/OneDrive - Technische Universität Ilmenau/A_UNI/BeachlorArbeit/Python_FRAP_Auswertung/UsedData/'

Loop = ['Control', 'CK666']
x = 0  # 'Control' data folder
x = 1 # 'CK666' data folder


"""____________________________________________________________________________
Importing Fitted Data (CSV) and Calculate Mean, Median & Standard Derivatian
read all files with extension ..."PhotofadeCorr.csv" or ..."noCorr.csv" """


filenames = glob.glob( Path + Loop[x] +'/Result Data/'+ "\*PhotofadeCorr.csv") #Photofade Corrected 
# filenames = glob.glob( Path + Loop[x] +'/Result Data/'+ "\*noCorr.csv")



""" for loop to iterate over all csv files and import them as list"""
fileList = []
for file in filenames:
    # reading csv files
    print(' ')
    print("Reading file = ",file)
    fileList.append(pd.read_csv(file)) # ad file to the file list
    
    
""" Extract data from list as numpy array"""
fileList = np.array(fileList)
fileList.shape

L = len(filenames)     # amount of files
D = fileList.shape[1]  # number of rows in imported data file: 31


""" Put all data (not parameter name) in one numpy array"""
dataArr = np.zeros(shape=(L,D)) 
for i in range(L):
    dataArr[i,:] = fileList[i,:,1] 
    
    

""" ----------------- Calculate Statistics --------------"""

Means = np.mean(dataArr, axis=0)
Medians = np.median(dataArr, axis=0)

Stds = np.zeros(shape=D) 
for i in range(D):
    Stds[i] = statistics.stdev(dataArr[:,i])    # Standard derivation for a random sample
    
# Stds2 = np.std(dataArr, axis=0)               # Standard derivation for a population
# print(np.vstack((Stds,Stds2)).transpose() )   # both compared


""" ------------- Configurate Result Table --------------"""

results = np.array([ np.zeros(D), np.zeros(D), Means, Stds, Medians, np.zeros(D), np.zeros(D)], dtype='object')

results[1,:] = (['','IE','I_immobile','y0','F1 ','A1','t1','t1_half','k1','F2','A2','t2','t2_half','k2','F_mobile','F_immobil',' ',' ','y0','A','t0',' ',' ','y0','A','t1','alpha',' ','m','x1','x2'])
results[:,16], results[:,21],results[0,:],results[5,:],results[6,:] = ['','','','','','',''],['','','','','','',''], ['']*D, ['']*D, ['']*D
results[0,1],results[0,18],results[0,23] = 'y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)','Opt.: y0 + A * (1 - exp(-t/τ))','Opt.: y0 + A * (1 - exp((-t/τ)^a))'
results[2,0],results[3,0],results[4,0] = 'Mean','Std','Median'
results[5,1],results[5,2] = 'No. of curves: ', L
results[6,0] = 'Used Curve Results:'


results = np.vstack((results , dataArr)) # Add all Data Arrays Unerneath the Result Array

results = np.delete(results, (17,22), 1) # delete empty columns



"""------------ Save Results as Excel Table ------------------------"""

dr = pd.DataFrame(data = results)  
dr.to_excel(Path + 'Result-'+ Loop[x]+'.xlsx', index=(None))



"""-------------------------------
Save (Seperatly) Mean Fit Parameters from  y_2 + A_1 *(1- np.exp(-t / t_1)) + A_2 *(1- np.exp(-t / t_2))
 y_2, A_1,t_1,A_2,t_2:"""
    
meanFitParams = np.array([ results[2,3] , results[2,5] , results[2,6] , results[2,10] , results[2,11] ])
np.savetxt(Path + Loop[x] + '/' + Loop[x] + '-Mean Params of 2-Exp-Fct.txt', X = meanFitParams, fmt ='%1.4f' ,delimiter=',')



###############################################################################

print('___________')
print('done')