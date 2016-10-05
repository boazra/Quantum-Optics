# -*- coding: utf-8 -*-
"""
Created on Wed Oct 05 19:49:11 2016

@author: Boaz
"""

from UDipole import Udipole
from numpy import min,argmin,linspace,diff,sqrt,pi
from matplotlib.pylab import figure,plot,legend,grid,axis,xlabel,ylabel,title,rc

h_bar = 1.054e-34
c3 = 4.9 #eV*(angstram)**3 http://www.physics.arizona.edu/~cronin/Research/vdw#20web/c3_theory.html
Kb = 1.38e-23

lambda_red = 1085 
lambda_blue = 729 
#c3 = (1.6e-19)*3.397e-35 #Serge
c3 =  1e6*(1.6e-19*c3*(1e-10)**3)*(1e9)**3/Kb # [microKelvin*nm**3]
lambda_res_1 = 780
lambda_res_2 = 795
input_power_red = 1000e-3 #[mW] need to be up to 40microW = 4e-2mW
input_power_blue = 330e-3 #[mW]
mode_area = 300e-8/(40*pi)  #[cm**2] 300e-8/(40*pi)

x = linspace(33,0.7*lambda_res_1,10000)
U_red_1 = Udipole(lambda_red,lambda_res_1,input_power_red,mode_area,x)
U_blue_1 = Udipole(lambda_blue,lambda_res_1,input_power_blue,0.67*mode_area,x)
U_red_2 = Udipole(lambda_red,lambda_res_2,input_power_red,mode_area,x)
U_blue_2 = Udipole(lambda_blue,lambda_res_2,input_power_blue,0.67*mode_area,x)
U_vdv = -c3/(x**3)

U_tot = U_red_1+U_blue_1+U_red_2+U_blue_2+U_vdv

dx = x[2]-x[1]
m = 1.42e-25
d2U_tot = diff(diff(1e6*Kb*U_tot)/dx)/dx
osc_freq = sqrt(abs(d2U_tot)/m)/(2*pi)
m = min(U_tot[100:])
i = argmin(U_tot[100:])

print osc_freq[i]

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':25})
rc('text', usetex=True)

figure(1)
plot(x,U_red_1,'r',x,U_blue_1,'b',x,U_red_2,'r',x,U_blue_2,'b',x,U_vdv,'g',x,U_tot,'k',linewidth=2.0)
legend(['red780','blue780','red795','blue795','vdv','tot'])
grid(True)

figure(2)
plot(x,U_red_1+U_red_2,'r',x,U_blue_1+U_blue_2,'b',x,U_vdv,'g',x,U_tot,'k',linewidth=2.0)
legend(['red','blue','vdv','tot'])
grid(True)
axis([0, 250, -1.1*max(U_blue_2+U_blue_1), 1.1*max(U_blue_2+U_blue_1)])

figure(3)
plot(x,U_tot,'k',linewidth=2.0)
xlabel('Distance from the Surface [nm]')
ylabel(r"Potential [$\mu$K]")
title('Dual-wavelength evanescent trap')
grid(True)
