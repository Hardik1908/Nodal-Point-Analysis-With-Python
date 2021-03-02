#Importing the libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

# Given input 
Pr = 6000  #reservoir pressure
D = 3.5    # tubing internal dia in inch
Choke_size = 64
PI = 1   #above bubble point
GLR = 1000
Water_cut = 0.25
Oil_API_gravity = 30
Water_specific_gravity = 1.05
Gas_specific_gravity= 0.65
Choke_constant= 10
GLR_exponent= 0.546
Choke_size_exponent = 1.89
Bw = 1
Wellhead_temperature = 100
Tubing_shoe_depth = 12000
Bottom_hole_temperature = 150
z = 0.9
air_density = 0.0765

# qo and pwf table
Pwf_ = Pr
qo = np.zeros((10,1))
Pwf = np.zeros((10,1))
for i in range (0,10):
  Pwf_  = Pwf_ - (6000/10)
  qo[i] = PI*(Pr - Pwf_)
  Pwf[i] = Pwf_

qo_Pwf = np.column_stack((qo,Pwf))

# Specific densities
Yo = 141.5/(131.5+ Oil_API_gravity)
Yg = Gas_specific_gravity
Yw = Water_specific_gravity

# iterating over Qo
k=0
PWH = np.zeros((10,1))
PWH_CPR = np.zeros((10,1))
for k,j in enumerate(qo_Pwf) :
  Qo = j[0]
  PWF = j[1]
  
  # rates 
  gas_rate = GLR*Qo
  water_rate = Water_cut * Qo
  oil_rate = Qo - water_rate
  
  #ratios
  WOR = water_rate/oil_rate
  GOR = gas_rate/oil_rate
  
  # M value
  M = 350.17*(Yo + WOR*Yw) + GOR*air_density*Yg
  
  
  # Iterating over assumed Pwh for given Qo and Pwf until  error is zero 
  error = 6
  dif = 0
  for i in range(1,100000):
    dif = dif + 0.1
    Pwh = PWF - dif
    # wellhead parameters
    Rs_head = Yg*(((Pwh/18)*(10**(0.0125*Oil_API_gravity)/10**(0.0009*Wellhead_temperature)))**1.2048)
    Bo_head = 0.971 + 0.000147*((Rs_head*((Yg/Yo)**0.5) + 1.25*Wellhead_temperature)**1.175)
    Vm_head = 5.615*(Bo_head + WOR*Bw) + (GOR - Rs_head)*(14.7/Pwh)*((Wellhead_temperature+460)/520)*(z/1)
    Rho1 = M/Vm_head
    
    
    # Bottom hole parameter
    Rs_bottom = Yg*(((PWF/18)*(10**(0.0125*Oil_API_gravity)/10**(0.0009*Bottom_hole_temperature)))**1.2048)
    Bo_bottom = 0.971 + 0.000147*((Rs_bottom*((Yg/Yo)**0.5) + 1.25*Bottom_hole_temperature)**1.175)
    Vm_bottom = 5.615*(Bo_bottom + WOR*Bw) + (GOR - Rs_bottom)*(14.7/PWF)*((Bottom_hole_temperature+460)/520)*(z/1)
    Rho2 = M/Vm_bottom
    
    # rho mixture
    Rho_mix = (Rho1 + Rho2)/2
    
    # Inertial force
    Inertial_force = 1.4737*(10**-5)*M*oil_rate/(D/12)
   
    # Friction factor
    Friction_factor = 4*10**(1.444-(2.5*math.log10(Inertial_force)))
   
    # Frition term
    K_bar = (Friction_factor*(oil_rate**2)*(M**2))/(7.4137*(10**10)*((D/12)**5))
   
    # error
    error = (144*(PWF-Pwh))/(Rho_mix+(K_bar/Rho_mix)) - Tubing_shoe_depth
    
    if -5  < error < 5:
      break
  # Calculating Pwh for CPR
  Pwh_CPR =  (Choke_constant*(GLR**GLR_exponent)*Qo)/(Choke_size**Choke_size_exponent)
  PWH_CPR[k] = Pwh_CPR 
  PWH_CPR = np.around(PWH_CPR)
  
  # making PWH vector
  PWH[k] = Pwh
  PWH[8] = 0  
  PWH = np.around(PWH)
  
print(PWH)
print(PWH_CPR)  
    
# Visualizing Result
final_table = np.column_stack((qo_Pwf, PWH, PWH_CPR))

plt.plot(final_table[[0,1,2,3,4,5,6,7,8],0], final_table[[0,1,2,3,4,5,6,7,8],2])
plt.plot(final_table[[0,1,2,3,4,5,6,7,8],0], final_table[[0,1,2,3,4,5,6,7,8],3]) 
plt.show()   
    
    
