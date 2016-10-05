# -*- coding: utf-8 -*-
"""
Created on Wed Oct 05 19:49:39 2016

@author: Boaz
"""

from numpy import pi,exp

def Udipole(lambda_detuned,lambda_res,input_power,mode_area,x):
    #   Calculates the dipole potential created by applying detuned laser on
    #   our cavity
    #   Units - lambdas in nm, input_power in mW, mode_area in cm**2
    #   x array is the distance from cavity in m
    
    Kb = 1.38e-23 #Joule/Kelvin
    finesse = 2.5e4
    Isat = 2.503 #Steck for cycling tansition
    
    gamma = 6e6        #[Hz] in the 1*gamma formulations (books, not ofer)
    c = 3e8            #[m/sec]
    h_bar = 1.054e-34  #[J*sec]
    
    res_freq = 2*pi*c/(lambda_res*1e-9)
    delta = 2*pi*c/(lambda_detuned*1e-9)-res_freq
    k = 2*pi/(lambda_detuned)
    
    I = (input_power/mode_area)*(0.3**2)*finesse*exp(-2*k*x)
    s = (I/Isat)*(1/(1+(2*delta/gamma)**2))
    
    U = h_bar*s*delta/(2*(1+s))
    U = 1e6*U/Kb # [microK]
    return U
