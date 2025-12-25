import numpy as np
import glob
from scipy.optimize import curve_fit
import os
import matplotlib
matplotlib.use('TkAgg')  # o prueba 'Agg' si no vas a mostrar ventanas
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
l = 1101  # número de filas a saltar en la lectura de los archivos CSV        
f = 3321 - l  
window_length = 151
polyorder = 3
############## todo lo que tiene terminacion _filt es una version filtrada de los datos
########################### DR G1 30G ###############################################
dr_g1_30g = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/DR-G1-30G.csv", skiprows=l, sep = ";", decimal = ",")
data_dr_g1_30g = dr_g1_30g.iloc[0:len(dr_g1_30g),0:2]
strlamda_dr_g1_30g = data_dr_g1_30g.iloc[:f,0].to_numpy()
strintensidad_dr_g1_30g = data_dr_g1_30g.iloc[:f,1].to_numpy()
lamda_dr_g1_30g = np.array([float(i) for i in strlamda_dr_g1_30g])
intensidad_dr_g1_30g = np.array([float(i) for i in strintensidad_dr_g1_30g])
intensidad_dr_g1_30g_filt = savgol_filter(intensidad_dr_g1_30g, window_length, polyorder)



############################# DR 30 ##########################################################
dr_30 = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/dr30.csv", skiprows=l, sep = ";", decimal = ",")
data_dr_30 = dr_30.iloc[0:len(dr_30),0:2]
strlamda_dr_30 = data_dr_30.iloc[:f,0].to_numpy()
strintensidad_dr_30 = data_dr_30.iloc[:f,1].to_numpy()
lamda_dr_30 = np.array([float(i) for i in strlamda_dr_30])
intensidad_dr_30 = np.array([float(i) for i in strintensidad_dr_30])
intensidad_dr_30_filt = savgol_filter(intensidad_dr_30, window_length, polyorder)



############################# R 30G ##########################################################
r_30g = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/R-30G.csv", skiprows=l, sep = ";", decimal = ",")
data_r_30g = r_30g.iloc[0:len(r_30g),0:2]
strlamda_r_30g = data_r_30g.iloc[:f,0].to_numpy()
strintensidad_r_30g = data_r_30g.iloc[:f,1].to_numpy()
lamda_r_30g = np.array([float(i) for i in strlamda_r_30g])
intensidad_r_30g = np.array([float(i) for i in strintensidad_r_30g])
intensidad_r_30g_filt = savgol_filter(intensidad_r_30g, window_length, polyorder)



############################# R 30 ##########################################################
r_30 = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/r30.csv", skiprows=l, sep = ";", decimal = ",")
data_r_30 = r_30.iloc[0:len(r_30),0:2]
strlamda_r_30 = data_r_30.iloc[:f,0].to_numpy()
strintensidad_r_30 = data_r_30.iloc[:f,1].to_numpy()
lamda_r_30 = np.array([float(i) for i in strlamda_r_30])
intensidad_r_30 = np.array([float(i) for i in strintensidad_r_30])
intensidad_r_30_filt = savgol_filter(intensidad_r_30, window_length, polyorder)

#R30 = (intensidad_r_30g - intensidad_dr_g1_30g)/(intensidad_r_30 - intensidad_dr_30)
R30 = (intensidad_r_30g_filt - intensidad_dr_g1_30g_filt)/(intensidad_r_30_filt - intensidad_dr_30_filt)

















############################## 45G ####################################
########################### DR G1 45G ###############################################
dr_g1_45g = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/DR-G1-45G.csv", skiprows=l, sep = ";", decimal = ",")
data_dr_g1_45g = dr_g1_45g.iloc[0:len(dr_g1_45g),0:2]
strlamda_dr_g1_45g = data_dr_g1_45g.iloc[:f,0].to_numpy()
strintensidad_dr_g1_45g = data_dr_g1_45g.iloc[:f,1].to_numpy()
lamda_dr_g1_45g = np.array([float(i) for i in strlamda_dr_g1_45g])
intensidad_dr_g1_45g = np.array([float(i) for i in strintensidad_dr_g1_45g])
intensidad_dr_g1_45g_filt = savgol_filter(intensidad_dr_g1_45g, window_length, polyorder)



############################# DR 45 ##########################################################
dr_45 = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/dr45.csv", skiprows=l, sep = ";", decimal = ",")
data_dr_45 = dr_45.iloc[0:len(dr_45),0:2]
strlamda_dr_45 = data_dr_45.iloc[:f,0].to_numpy()
strintensidad_dr_45 = data_dr_45.iloc[:f,1].to_numpy()
lamda_dr_45 = np.array([float(i) for i in strlamda_dr_45])
intensidad_dr_45 = np.array([float(i) for i in strintensidad_dr_45])
intensidad_dr_45_filt = savgol_filter(intensidad_dr_45, window_length, polyorder)



############################# R 45G ##########################################################
r_45g = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/R-45G.csv", skiprows=l, sep = ";", decimal = ",")
data_r_45g = r_45g.iloc[0:len(r_45g),0:2]
strlamda_r_45g = data_r_45g.iloc[:f,0].to_numpy()
strintensidad_r_45g = data_r_45g.iloc[:f,1].to_numpy()
lamda_r_45g = np.array([float(i) for i in strlamda_r_45g])
intensidad_r_45g = np.array([float(i) for i in strintensidad_r_45g])
intensidad_r_45g_filt = savgol_filter(intensidad_r_45g, window_length, polyorder)



############################# R 45 ##########################################################
r_45 = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/r45.csv", skiprows=l, sep = ";", decimal = ",")
data_r_45 = r_45.iloc[0:len(r_45),0:2]
strlamda_r_45 = data_r_45.iloc[:f,0].to_numpy()
strintensidad_r_45 = data_r_45.iloc[:f,1].to_numpy()
lamda_r_45 = np.array([float(i) for i in strlamda_r_45])
intensidad_r_45 = np.array([float(i) for i in strintensidad_r_45])
intensidad_r_45_filt = savgol_filter(intensidad_r_45, window_length, polyorder)

#R45 = (intensidad_r_45g - intensidad_dr_g1_45g)/(intensidad_r_45 - intensidad_dr_45)
R45 = (intensidad_r_45g_filt - intensidad_dr_g1_45g_filt)/(intensidad_r_45_filt - intensidad_dr_45_filt)



############################## 60G #################################
########################### DR G1 60G ###############################################
dr_g1_60g = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/DR-G1-60G.csv", skiprows=l, sep = ";", decimal = ",")
data_dr_g1_60g = dr_g1_60g.iloc[0:len(dr_g1_60g),0:2]
strlamda_dr_g1_60g = data_dr_g1_60g.iloc[:f,0].to_numpy()
strintensidad_dr_g1_60g = data_dr_g1_60g.iloc[:f,1].to_numpy()
lamda_dr_g1_60g = np.array([float(i) for i in strlamda_dr_g1_60g])
intensidad_dr_g1_60g = np.array([float(i) for i in strintensidad_dr_g1_60g])
intensidad_dr_g1_60g_filt = savgol_filter(intensidad_dr_g1_60g, window_length, polyorder)



############################# DR 60 ##########################################################
dr_60 = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/dr60.csv", skiprows=l, sep = ";", decimal = ",")
data_dr_60 = dr_60.iloc[0:len(dr_60),0:2]
strlamda_dr_60 = data_dr_60.iloc[:f,0].to_numpy()
strintensidad_dr_60 = data_dr_60.iloc[:f,1].to_numpy()
lamda_dr_60 = np.array([float(i) for i in strlamda_dr_60])
intensidad_dr_60 = np.array([float(i) for i in strintensidad_dr_60])
intensidad_dr_60_filt = savgol_filter(intensidad_dr_60, window_length, polyorder)



############################# R 60G ##########################################################
r_60g = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/R-60G.csv", skiprows=l, sep = ";", decimal = ",")
data_r_60g = r_60g.iloc[0:len(r_60g),0:2]
strlamda_r_60g = data_r_60g.iloc[:f,0].to_numpy()
strintensidad_r_60g = data_r_60g.iloc[:f,1].to_numpy()
lamda_r_60g = np.array([float(i) for i in strlamda_r_60g])
intensidad_r_60g = np.array([float(i) for i in strintensidad_r_60g])
intensidad_r_60g_filt = savgol_filter(intensidad_r_60g, window_length, polyorder)



############################# R 60 ##########################################################
r_60 = pd.read_csv("C:/Users/Asus/caracterizacion_pelicula_delgada/r60.csv", skiprows=l, sep = ";", decimal = ",")
data_r_60 = r_60.iloc[0:len(r_60),0:2]
strlamda_r_60 = data_r_60.iloc[:f,0].to_numpy()
strintensidad_r_60 = data_r_60.iloc[:f,1].to_numpy()
lamda_r_60 = np.array([float(i) for i in strlamda_r_60])
intensidad_r_60 = np.array([float(i) for i in strintensidad_r_60])
intensidad_r_60_filt = savgol_filter(intensidad_r_60, window_length, polyorder)

#R60 = (intensidad_r_60g - intensidad_dr_g1_60g)/(intensidad_r_60 - intensidad_dr_60)
R60 = (intensidad_r_60g_filt - intensidad_dr_g1_60g_filt)/(intensidad_r_60_filt - intensidad_dr_60_filt)












reflectance_data = np.column_stack((lamda_r_30, R30, R45, R60))
np.savetxt("C:/Users/Asus/caracterizacion_pelicula_delgada/filt_reflectance_data.csv", reflectance_data, delimiter = ",", header='wavelength_nm,R30,R45,R60', comments='')