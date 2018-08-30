"""
This is to calculate GHG species RF based on the IPCC formulate
"""

import numpy as np
import math
def CO2_RF(C2,C1):
	##c2 and c1 are  CO2 mixing ratio in ppm
	return 5.35*np.log(C2*1.0/C1)

def CFC11_RF(C2,C1):
	##c2 and c1 are CFC11 mixing ratio in ppb
	return 0.25*(C2-C1)
	
def CFC12_RF(C2,C1):
	##c2 and c1 are CFC12 mixing ratio in ppb
	return 0.32*(C2-C1)
	
def f(M,N):
	#### A FUNCITON FOR CH4 AND N2O RF CALCULATION
	#### M,N are CH4 and N2O in ppb
	return 0.47*np.log(1+2.01*10**(-5)*np.power(M*N,0.75)+5.31*10**(-51)*M*np.power(M*N,1.52))
	
def CH4_RF(CH4_2,CH4_1,N2O):
	##c2 and c1 are CH4 mixing ratio in ppb
	return 0.036*(math.sqrt(CH4_2)-math.sqrt(CH4_1))-(f(CH4_2,N2O)-f(CH4_1,N2O))

def N2O_RF(N2O_2,N2O_1,CH4):
	##c2 and c1 are N2O mixing ratio in ppb
	return 0.12*(math.sqrt(N2O_2)-math.sqrt(N2O_1))-(f(CH4,N2O_2)-f(CH4,N2O_1))
	
def calculate_print_RF():
	### These values corresponde to 1970 and 2010
	CO2_1 = 324.985;   CO2_2 = 389.32416;   #PPM
	CH4_1 = 1385.75;   CH4_2 = 1778.6749;   #PPB
	N2O_1 = 295.2 ;    N2O_2 = 323.06113;   #PPB
	CFC11_1 = 159.24573/1000; CFC11_2 = 747.47797/1000;  #PPB
	CFC12_1 = 107.45/1000; CFC12_2 = 524.93396/1000;  #PPB
	
	CO2= CO2_RF(CO2_2,CO2_1)
	CH4= CH4_RF(CH4_2,CH4_1,N2O_2)
	N2O= N2O_RF(N2O_2,N2O_1,CH4_2)
	CFC11= CFC11_RF(CFC11_2,CFC11_1)
	CFC12= CFC12_RF(CFC12_2,CFC12_1)
	print 'C02',CO2,'CH4',CH4,'N2O',N2O,'CFC11',CFC11,'CFC12',CFC12
	print 'GHG', CO2+CH4+N2O+CFC11+CFC12
	
	
	
calculate_print_RF()
	
	
	
	
	
	
	
	