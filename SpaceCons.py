""" Constants in Aerospace manner """

import numpy as np

# area to mass (A/m) of spacecraft, m^2/kg
area2mass = 0.1 # m^2/kg
REFLECTABILITY = 0.8

# uEarth, m^3/s; u = GM
MU_EARTH = 398600.4415e9 

# gravitational constant, m^3/kg*s^2
G = 6.673e-11 
# Earth mass: kg 
MASS_EARTH = 5.973332e24 
# mean equatorial radius of Earth, JGM-3 model, m 
RADIUS_EARTH = 6378.1363e3 
# rotational velocity of Earth, rad/s
OMEGA_EARTH = 7.292115e-5 

# GEO height, m
HGEO = 35786.7e3
# GEO radius, m
RGEO = HGEO + RADIUS_EARTH

### Earth gravitational coefficients
# coefficients, have no unit
C00 =  1.0
C10 =  0.0
C20 = -1.083e-3 
C21 =  0.0	   
C22 =  1.574e-6
C30 =  2.532e-6 
C31 =  2.192e-6  
C32 =  3.090e-7  
C33 =  1.005e-7 
C40 =  1.620e-6 
C41 = -5.088e-7 
C42 =  7.842e-8  
C43 =  5.921e-8  
C44 = -3.984e-9 

J2 = -C20 
J3 = -C30  
J4 = -C40

S22 = -9.038e-7
S31 =  2.684e-7  
S32 = -2.114e-7  
S33 =  1.972e-7 
S41 = -4.491e-7  
S42 =  1.482e-7  
S43 = -1.201e-8  
S44 =  6.526e-9

### for Sun, P1043 Vallado, V4
MU_SUN	   = 1.32712428e20  # m^3/s
AU		   = 1.4959787e11   # m
RADIUS_SUN = 6.96e8		 # m

### for Moon, P1041 Vallado
MU_MOON	   = 4902.799e9  # m^3/s
Earth_Moon = 384400e3	# m

### for solar radiation pressure in Earth orbit
P_SUN = 4.56e-6   # N/m^2


### for realistic sail
# solar sail book P50
# Square sail & Heliogyro, P50, the Solarsail Book
REALSAIL_R_TILD  = 0.88
REALSAIL_S_TILD  = 0.94
REALSAIL_EPSTONF = 0.05
REALSAIL_EPSTONB = 0.55
REALSAIL_BF	     = 0.79
REALSAIL_BB	     = 0.55

### coordinates of ECI
G1 = np.array([1, 0, 0])
G2 = np.array([0, 1, 0])
G3 = np.array([0, 0, 1])

### time constants
TIMEZONE = 0.0
# July 1st 2017
DUT1 = 0.359485
DAT  = 37.0
# # July 1st 2012
# DUT1 = 0.413223
# DAT  = 35.0
# # July 1st 2007
# DUT1 = -0.157336
# DAT  = 33.0











