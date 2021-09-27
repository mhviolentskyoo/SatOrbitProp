"""
Orbital disturbance functions, includes:
class ParDeri		:  includes partial derivatives
						U2acc, partial derivatives of UJ2 UJ22 UJ3 UJ31 UJ32 UJ33
def acc(posECEF)	:  return J2 J22 J3 J31 J32 J33 accelerations in ECEF

"""

import numpy as np
import pdb
from scipy.special import eval_legendre

import SpaceCons



""" Partial derivatives """
class ParDeri:
	def __init__(self):
		pass

	@staticmethod
	def par_U2acc(posECEF, parU):
		"""
			Use partial direvatives of gravitational potential to calculate acceleration
			Ref   : 
					Vallado, V4, P550, Eq.8-27, don't consider two-body acceleration
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda
			outputs: 
					accECEF(1x3 numpy array): [ai,aj,ak]
						acceleration in ECEF
		"""
		r = np.linalg.norm(posECEF)
		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]
		parU_r, parU_phi, parU_lam = parU[0], parU[1], parU[2]

		# print("radius = ", r)
		# print("posECEF = ", ri, rj ,rk)
		# print("parU = ", parU)
		# pdb.set_trace() # break ponit

		ai = (1/r * parU_r - rk / r**2 / (ri**2+rj**2)**(1/2) * parU_phi) * ri \
			 - (1/(ri**2+rj**2) * parU_lam) * rj
		aj = (1/r * parU_r - rk / r**2 / (ri**2+rj**2)**(1/2) * parU_phi) * rj \
			 + (1/(ri**2+rj**2) * parU_lam) * ri
		ak = 1/r * parU_r * rk + (ri**2+rj**2)**(1/2) / r**2 * parU_phi

		accECEF = np.array([ai, aj, ak])

		# print("accECEF = ", accECEF)
		# pdb.set_trace() # break ponit

		return accECEF

	@staticmethod
	def par_UJ2(posECEF):
		"""
			Use J2 potential to calculate partial derivatives
			Ref   : 
					Vallado, V4, P550
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			output: 
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda
		"""
		r = np.linalg.norm(posECEF)

		# print(posECEF, r)
		# pdb.set_trace() # break ponit

		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]

		parU_r   = 3/2 * SpaceCons.MU_EARTH * SpaceCons.J2 \
				   * (SpaceCons.RADIUS_EARTH/r**2)**2 * (3.0 * rk**2/r**2 - 1)
		parU_phi = -3.0 * SpaceCons.MU_EARTH * SpaceCons.J2 \
				   * (SpaceCons.RADIUS_EARTH)**2 * 1/r**5 * rk * (ri**2+rj**2)**(1/2)
		parU_lam = 0

		parU = np.array([parU_r, parU_phi, parU_lam])

		# print("parUJ2 = ", parU)
		# pdb.set_trace() # break ponit

		return parU

	@staticmethod
	def par_UJ22(posECEF):
		"""
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			output: 
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda

		"""
		r = np.linalg.norm(posECEF)
		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]

		parU_r   = -9.0 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**2 / r**6 \
					* (SpaceCons.C22*(ri**2-rj**2) + SpaceCons.S22*2*ri*rj)
		parU_phi = -6.0 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**2 / r**5 \
					* rk * (SpaceCons.C22*(ri**2-rj**2) + SpaceCons.S22*2*ri*rj) / (ri**2+rj**2)**(1/2)
		parU_lam = 3.0 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**2 / r**5 \
					* (-4*SpaceCons.C22*ri*rj + 2*SpaceCons.S22*(ri**2-rj**2))

		parU = np.array([parU_r, parU_phi, parU_lam])

		# print("parUJ22 = ", parU)
		# pdb.set_trace() # break ponit

		return parU

	@staticmethod
	def par_UJ3(posECEF):
		"""
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			output: 
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda
		"""
		r = np.linalg.norm(posECEF)
		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]

		coslam = ri / (ri**2 + rj**2)**(1/2)
		sinlam = rj / (ri**2 + rj**2)**(1/2)
		cosphi = (ri**2 + rj**2)**(1/2) / r
		sinphi = rk / r

		parU_r   = 2 * SpaceCons.MU_EARTH * SpaceCons.J3 * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**5 * (5 * sinphi**3 - 3 * sinphi)
		parU_phi = -0.5 * SpaceCons.MU_EARTH * SpaceCons.J3 * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * (15 * sinphi**2 * cosphi - 3 * cosphi)
		parU_lam = 0

		parU = np.array([parU_r, parU_phi, parU_lam])

		# print("parUJ3 = ", parU)
		# pdb.set_trace() # break ponit

		return parU

	@staticmethod
	def par_UJ31(posECEF):
		"""
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			output: 
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda
		"""
		r = np.linalg.norm(posECEF)
		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]

		coslam = ri / (ri**2 + rj**2)**(1/2)
		sinlam = rj / (ri**2 + rj**2)**(1/2)
		cosphi = (ri**2 + rj**2)**(1/2) / r
		sinphi = rk / r

		parU_r   = -2 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**5 * cosphi * (15 * sinphi**2 - 3) \
				   * (SpaceCons.C31 * coslam + SpaceCons.S31 * sinlam)
		parU_phi = 0.5 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * (SpaceCons.C31 * coslam + SpaceCons.S31 * sinlam) \
				   * (30 * sinphi * cosphi**2 - 15 * sinphi**3 + 3 * sinphi)
		parU_lam = 0.5 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * cosphi * (15 * sinphi**2 - 3) \
				   * (-SpaceCons.C31 * sinphi + SpaceCons.S31 * cosphi)

		parU = np.array([parU_r, parU_phi, parU_lam])

		# print("parUJ31 = ", parU)
		# pdb.set_trace() # break ponit

		return parU

	@staticmethod
	def par_UJ32(posECEF):
		"""
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			output: 
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda
		"""
		r = np.linalg.norm(posECEF) 
		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]

		coslam  = ri / (ri**2 + rj**2)**(1/2)
		sinlam  = rj / (ri**2 + rj**2)**(1/2)
		cosphi  = (ri**2 + rj**2)**(1/2) / r
		sinphi  = rk / r
		cos2lam = coslam**2 - sinlam**2
		sin2lam = 2 * sinlam * coslam

		parU_r   = -4 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**5 * 15 * cosphi**2 * sinphi \
				   * (SpaceCons.C32 * cos2lam + SpaceCons.S32 * sin2lam)
		parU_phi = 15 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * (SpaceCons.C32 * cos2lam + SpaceCons.S32 * sin2lam) \
				   * (-2 * cosphi * sinphi**2 + cosphi**3)
		parU_lam = 15 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * cosphi**2 * sinphi \
				   * (-2 * SpaceCons.C32 * sin2lam + 2 * SpaceCons.S32 * cos2lam)

		parU = np.array([parU_r, parU_phi, parU_lam])

		# print("parUJ32 = ", parU)
		# pdb.set_trace() # break ponit

		return parU

	@staticmethod
	def par_UJ33(posECEF):
		"""
			inputs: 
					posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			output: 
					parU(1x3 numpy array): [parU_r, parU_phi, parU_lam]
						partial direvatives of potential U to r,phi,lambda
		"""
		r = np.linalg.norm(posECEF) 
		ri, rj, rk = posECEF[0], posECEF[1], posECEF[2]

		coslam  = ri / (ri**2 + rj**2)**(1/2)
		sinlam  = rj / (ri**2 + rj**2)**(1/2)
		cosphi  = (ri**2 + rj**2)**(1/2) / r
		sinphi  = rk / r
		cos3lam = 4 * coslam**3 - 3 * coslam
		sin3lam = 3 * sinlam - 4 * sinlam**3

		parU_r   = -4 * SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**5 * 15 * cosphi**3 \
				   * (SpaceCons.C33 * cos3lam + SpaceCons.S33 * sin3lam)
		parU_phi = SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * 15 * 3 * cosphi**2 * (-sinphi) \
				   * (SpaceCons.C33 * cos3lam + SpaceCons.S33 * sin3lam)
		parU_lam = SpaceCons.MU_EARTH * (SpaceCons.RADIUS_EARTH)**3 \
				   / r**4 * 15 * cosphi**3 \
				   * (-3 * SpaceCons.C33 * sin3lam + 3 * SpaceCons.S33 * cos3lam)

		parU = np.array([parU_r, parU_phi, parU_lam])

		# print("parUJ33 = ", parU)
		# pdb.set_trace() # break ponit

		return parU



""" J2 acceleration """
def J2_pert(posECEF):
	"""
	Calculate J2 perturbation using position in ECEF
	inputs:
			posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			in m
	outputs:
			J2_acc(1x3 numpy array): J2 acceleration in ECEF [ai, aj, ak]
			in m/s
	"""
	# calculate partial derivatives for UJ2
	parU_J2 = ParDeri.par_UJ2(posECEF)

	# find the J2 acceleration in ECEF
	J2_acc  = ParDeri.par_U2acc(posECEF, parU_J2)

	# print("J2_acc = ", J2_acc)
	# pdb.set_trace() # break ponit

	return J2_acc



""" J22 acceleration """
def J22_pert(posECEF):
	"""
	Calculate J22 perturbation using position in ECEF
	inputs:
			posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			in m
	outputs:
			J22_acc(1x3 numpy array): J22 acceleration in ECEF [ai, aj, ak]
			in m/s
	"""
	# calculate partial derivatives for UJ22
	parU_J22 = ParDeri.par_UJ22(posECEF)

	# find the J2 acceleration in ECEF
	J22_acc = ParDeri.par_U2acc(posECEF, parU_J22)

	# print("J22_acc = ", J22_acc)
	# pdb.set_trace() # break ponit

	return J22_acc



""" J3 acceleration """
def J3_pert(posECEF):
	"""
	Calculate J3 perturbation using position in ECEF
	inputs:
			posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			in m
	outputs:
			J3_acc(1x3 numpy array): J3 acceleration in ECEF [ai, aj, ak]
			in m/s
	"""
	# calculate partial derivatives for UJ3
	parU_J3 = ParDeri.par_UJ3(posECEF)

	# find the J3 acceleration in ECEF
	J3_acc = ParDeri.par_U2acc(posECEF, parU_J3)

	# print("J3_acc = ", J3_acc)
	# pdb.set_trace() # break ponit

	return J3_acc



""" J31 acceleration """
def J31_pert(posECEF):
	"""
	Calculate J31 perturbation using position in ECEF
	inputs:
			posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			in m
	outputs:
			J31_acc(1x3 numpy array): J31 acceleration in ECEF [ai, aj, ak]
			in m/s
	"""
	# calculate partial derivatives for UJ31
	parU_J31 = ParDeri.par_UJ31(posECEF)

	# find the J31 acceleration in ECEF
	J31_acc = ParDeri.par_U2acc(posECEF, parU_J31)

	# print("J31_acc = ", J31_acc)
	# pdb.set_trace() # break ponit

	return J31_acc



""" J32 acceleration """
def J32_pert(posECEF):
	"""
	Calculate J32 perturbation using position in ECEF
	inputs:
			posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			in m
	outputs:
			J32_acc(1x3 numpy array): J32 acceleration in ECEF [ai, aj, ak]
			in m/s
	"""
	# calculate partial derivatives for UJ32
	parU_J32 = ParDeri.par_UJ32(posECEF)

	# find the J32 acceleration in ECEF
	J32_acc = ParDeri.par_U2acc(posECEF, parU_J32)

	# print("J32_acc = ", J32_acc)
	# pdb.set_trace() # break ponit

	return J32_acc



""" J33 acceleration """
def J33_pert(posECEF):
	"""
	Calculate J33 perturbation using position in ECEF
	inputs:
			posECEF(1x3 numpy array): coordinates in ECEF [ri,rj,rk]
			in m
	outputs:
			J33_acc(1x3 numpy array): J33 acceleration in ECEF [ai, aj, ak]
			in m/s
	"""
	# calculate partial derivatives for UJ33
	parU_J33 = ParDeri.par_UJ33(posECEF)

	# find the J33 acceleration in ECEF
	J33_acc = ParDeri.par_U2acc(posECEF, parU_J33)

	# print("J33_acc = ", J33_acc)
	# pdb.set_trace() # break ponit

	return J33_acc



""" Sun & Moon position vectors in ECI """
class thirdBodyPosECI:

	def __init__(self):
		pass

	@staticmethod
	def SunPosECI(tut1, ttdb):
		"""
		Calculate the Sun position vector in ECI
			inputs:
					tut1, ttdb
			outputs:
					Sun_posECI(1x3 numpy array): Sun position vector in ECI
					in m
			Ref: Vallado, V4, P279 & Vallado functions sun
				 make a subtle improvement of ttdb ~= tut1
		"""
		twopi	 = 2.0 * np.pi
		deg2rad  = np.pi / 180.0
		meanlong = 280.460 + 36000.77 * tut1
		meanlong = np.remainder(meanlong, 360.0)  # deg

		# Vallado: ttdb = tut1
		# we input ttdb directly here
		# note that the dif between ttdb & tut1 ~= 1e-8
		meananomaly = 357.5277233  + 35999.05034 * ttdb
		meananomaly = np.remainder(meananomaly * deg2rad, twopi)  # rad
		if meananomaly < 0.0 :
		   meananomaly = twopi + meananomaly

		eclplong = meanlong + 1.914666471 * np.sin(meananomaly) \
				  + 0.019994643 * np.sin(2.0 * meananomaly)  # deg
		eclplong = np.remainder(eclplong, 360.0)  # deg

		obliquity = 23.439291 - 0.0130042 * ttdb  # deg

		eclplong  = eclplong * deg2rad
		obliquity = obliquity * deg2rad

		magr = 1.000140612  - 0.016708617 * np.cos(meananomaly) \
			   - 0.000139589 * np.cos(2.0 * meananomaly)
	 
		ri = magr * np.cos(eclplong)
		rj = magr * np.cos(obliquity) * np.sin(eclplong)
		rk = magr * np.sin(obliquity) * np.sin(eclplong)

		# convert AU to m
		Sun_posECI = np.array([ri, rj, rk]) * SpaceCons.AU

		return Sun_posECI

	@staticmethod
	def MoonPosECI(ttdb):
		"""
		Calculate the Moon position vector in ECI
			inputs:
					ttdb
			outputs:
					Moon_posECI(1x3 numpy array): Moon position vector in ECI
					in m
			Ref: Vallado, V4, P279 & Vallado functions moon
		"""
		twopi	 = 2.0 * np.pi
		deg2rad  = np.pi / 180.0

		eclplong = 218.32  + 481267.8813 * ttdb \
				   + 6.29 * np.sin((134.9+477198.85*ttdb)*deg2rad) \
				   - 1.27 * np.sin((259.2-413335.38*ttdb)*deg2rad) \
				   + 0.66 * np.sin((235.7+890534.23*ttdb)*deg2rad) \
				   + 0.21 * np.sin((269.9+954397.70*ttdb)*deg2rad) \
				   - 0.19 * np.sin((357.5+35999.05*ttdb)*deg2rad) \
				   - 0.11 * np.sin((186.6+966404.05*ttdb)*deg2rad) # deg

		eclplat = 5.13 * np.sin( (93.3+483202.03*ttdb)*deg2rad ) \
				  + 0.28 * np.sin( (228.2+960400.87*ttdb)*deg2rad ) \
				  - 0.28 * np.sin( (318.3+6003.18*ttdb)*deg2rad ) \
				  - 0.17 * np.sin( (217.6-407332.20*ttdb)*deg2rad )  # deg

		hzparal = 0.9508 + 0.0518 * np.cos( (134.9+477198.85*ttdb)*deg2rad ) \
				  + 0.0095 * np.cos( (259.2-413335.38*ttdb)*deg2rad ) \
				  + 0.0078 * np.cos( (235.7+890534.23*ttdb)*deg2rad ) \
				  + 0.0028 * np.cos( (269.9+954397.70*ttdb)*deg2rad ) # deg

		eclplong = np.remainder( eclplong*deg2rad, twopi )
		eclplat  = np.remainder( eclplat*deg2rad, twopi )
		hzparal  = np.remainder( hzparal*deg2rad, twopi )

		obliquity = 23.439291  - 0.0130042 * ttdb # deg
		obliquity = obliquity * deg2rad

		# find the geocentric direction cosines
		l = np.cos(eclplat) * np.cos(eclplong)
		m = np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong) \
			- np.sin(obliquity) * np.sin(eclplat)
		n = np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong) \
			+ np.cos(obliquity) * np.sin(eclplat)

		# calculate moon position vector
		# small correction from Vallado: modify 1 to RE (in m)
		magr = SpaceCons.RADIUS_EARTH / np.sin(hzparal)
		ri = magr * l 
		rj = magr * m
		rk = magr * n

		Moon_posECI = np.array([ri, rj, rk])

		return Moon_posECI



""" thirdbody gravitational accelerations """
class thirdBodyAccECI:

	def __init__(self):
		pass

	@staticmethod
	def thirdbody_Sun(satPosECI, SunPosECI):
		"""
		Calculate the thirdbody (Sun) gravitational acceleration vector in ECI
			inputs:
					satPosECI(1x3 numpy array): satellite position vector in ECI
					in m
					SunPosECI(1x3 numpy array): Sun position vector in ECI
					in m
			outputs:
					acc_GraSun(1x3 numpy array): thirdbody (Sun) gravitational acceleration vector in ECI
					in m/s
			Ref: Vallado, V4, P574 
				 Mei H, Damaren C J, Zhan X. 
				 Hybrid removal of end-of-life geosynchronous satellites using solar radiation pressure and impulsive thrusts[J]. 
				 Advances in Space Research, 2020, 66(4): 974-991.
		"""
		# distance vector
		r_sat_Sun = SunPosECI - satPosECI
		# distances
		rEarthSun = np.linalg.norm(SunPosECI)
		rEarthSat = np.linalg.norm(satPosECI)
		cosinEpsilon = np.matmul(SunPosECI, satPosECI) / rEarthSun / rEarthSat
 
		# vector A (see the ref paper)
		vectorA = - SpaceCons.MU_SUN * satPosECI / rEarthSun**3
 
		# vector B (see the ref paper)
		cosinEpsilon = np.real(cosinEpsilon)
		tempB1 = eval_legendre(1, cosinEpsilon)
		B1 = tempB1 * (rEarthSat/rEarthSun)
		tempB2 = eval_legendre(2, cosinEpsilon)
		B2 = tempB2 * (rEarthSat/rEarthSun)**2
		# B terms and vector B 
		B = B1 + B2		 
		vectorB = SpaceCons.MU_SUN * r_sat_Sun / rEarthSun**3 \
				  * (3.0*B + 3.0*B**2 + B**3)
 
		# total
		acc_GraSun = vectorA + vectorB

		return acc_GraSun

	@staticmethod
	def thirdbody_Moon(satPosECI, MoonPosECI):
		"""
		Calculate the thirdbody (Moon) gravitational acceleration vector in ECI
			inputs:
					satPosECI(1x3 numpy array) : satellite position vector in ECI
					in m
					MoonPosECI(1x3 numpy array): Moon position vector in ECI
					in m
			outputs:
					acc_GraMoon(1x3 numpy array): thirdbody (Moon) gravitational acceleration vector in ECI
					in m/s
			Ref: Vallado, V4, P574 
				 Mei H, Damaren C J, Zhan X. 
				 Hybrid removal of end-of-life geosynchronous satellites using solar radiation pressure and impulsive thrusts[J]. 
				 Advances in Space Research, 2020, 66(4): 974-991.
		"""
		# distance vector
		r_sat_Moon = MoonPosECI - satPosECI
		# distances
		rEarthMoon = np.linalg.norm(MoonPosECI)
		rEarthSat  = np.linalg.norm(satPosECI)
		cosinEpsilon = np.matmul(MoonPosECI, satPosECI) / rEarthMoon / rEarthSat

		# vector A (see the ref paper)
		vectorA = - SpaceCons.MU_MOON * satPosECI / rEarthMoon**3
 
		# vector B (see the ref paper)
		cosinEpsilon = np.real(cosinEpsilon)
		tempB1 = eval_legendre(1, cosinEpsilon)
		B1 = tempB1 * (rEarthSat/rEarthMoon)
		tempB2 = eval_legendre(2, cosinEpsilon)
		B2 = tempB2 * (rEarthSat/rEarthMoon)**2
		tempB3 = eval_legendre(3, cosinEpsilon)
		B3 = tempB3 * (rEarthSat/rEarthMoon)**3
		# B terms and vector B 
		B = B1 + B2 + B3
		vectorB = SpaceCons.MU_MOON * r_sat_Moon / rEarthMoon**3 \
				  * (3.0*B + 3.0*B**2 + B**3)

		# total
		acc_GraMoon = vectorA + vectorB

		return acc_GraMoon




""" SRP related """
def controlAngles(xECI, Tut1, Ttdb):
	"""
	Calculate solar sail control angles
		inputs:
				xECI(1x6 numpy array) : satellite position & velocity vector in ECI [rx, ry, rz, vx, vy, vz]
				in m
				Tut1, Ttdb(float): 
		outputs:
				controlAngles(1x2 numpy array): sloar sail control angles [alpha, delta]
				in rad
				when sail moves towards the Sun: sail normal vector aligns with Sunline vector:
				SRP maximized
				when sail moves away from the Sun: sail normal vector perpendicular to Sunline vector:
				SRP = 0
		Ref: Vallado, V4, P574 
			 Mei H, Damaren C J, Zhan X. 
			 Hybrid removal of end-of-life geosynchronous satellites using solar radiation pressure and impulsive thrusts[J]. 
			 Advances in Space Research, 2020, 66(4): 974-991.
	"""
	# sun position vector in ECI
	Sun_posECI  = thirdBodyPosECI.SunPosECI(Tut1, Ttdb)

	# sunline vector in ECI
	r_sat_Sun = Sun_posECI - xECI[0:3]
	# opposite direction of r_sat_Sun, and normalize
	vec_u = - r_sat_Sun / np.linalg.norm(r_sat_Sun) 

	# nominal control angles, in rad
	alpha, delta = 0, np.pi # maximize SRP
	if np.matmul(vec_u, xECI[3:6]) <= 0:
		alpha, delta = np.pi/2, np.pi # SRP = 0

	sailControlAngles = np.array([alpha, delta])

	# print(sailControlAngles)
	# pdb.set_trace() # break ponit

	return sailControlAngles


def SRP_perfect(xECI, sailControlAngles, Tut1, Ttdb):
	"""
	Calculate SRP due to perfect solar sail in ECI
		inputs:
				xECI(1x6 numpy array) : satellite position & velocity vector in ECI [rx, ry, rz, vx, vy, vz]
				in m
				Tut1, Ttdb(float): 
		outputs:
				SRP_perfect(1x3 numpy array): SRP in ECI [ai, aj, ak]
				eclipse_ratio(float): eclipse condition, 1:light, 0:umbra, (0,1): punumbra
		Ref: Vallado, V4, P574 
			 Mei H, Damaren C J, Zhan X. 
			 Hybrid removal of end-of-life geosynchronous satellites using solar radiation pressure and impulsive thrusts[J]. 
			 Advances in Space Research, 2020, 66(4): 974-991.
	"""
	# sun position vector in ECI
	Sun_posECI  = thirdBodyPosECI.SunPosECI(Tut1, Ttdb)

	# sunline vector in ECI
	r_sat_Sun = Sun_posECI - xECI[0:3]
	# opposite direction of r_sat_Sun, and normalize
	vec_u = - r_sat_Sun / np.linalg.norm(r_sat_Sun) 


	### sail orientation (control)
	# Construct frame Fs (see ref paper)
	S1 = vec_u
	S3 = (SpaceCons.G3-np.matmul(SpaceCons.G3,vec_u)*vec_u) / np.linalg.norm(SpaceCons.G3-np.matmul(SpaceCons.G3,vec_u)*vec_u)
	S2 = np.cross(S3,S1) / np.linalg.norm(np.cross(S3,S1))
	# Cgs : rotation matrix from Fs to ECI
	Cgs = np.array([[np.matmul(SpaceCons.G1,S1), np.matmul(SpaceCons.G1,S2), np.matmul(SpaceCons.G1,S3)],
					[np.matmul(SpaceCons.G2,S1), np.matmul(SpaceCons.G2,S2), np.matmul(SpaceCons.G2,S3)],
					[np.matmul(SpaceCons.G3,S1), np.matmul(SpaceCons.G3,S2), np.matmul(SpaceCons.G3,S3)]])
	# coordinates of sail normal vector in Fs
	alpha, delta = sailControlAngles[0], sailControlAngles[1]
	ns = np.array([np.cos(alpha),
				   np.sin(alpha) * np.sin(delta),
				   np.sin(alpha) * np.cos(delta)])
	# coordinates of sail normal vector in ECI
	vec_n = np.matmul(Cgs, ns)


	### SRP
	cosinAlpha = np.matmul(vec_n, vec_u)
	tempP_sun = SpaceCons.P_SUN
	## eclipse (see ref)
	# eclipse_ratio: 1:light, 0:umbra, (0,1): punumbra
	eclipse_ratio = 1
	if np.matmul(Sun_posECI, xECI[0:3]) < 0: # when the sat is behind the Earth
		# distances: Earth2Sun, Earth2Sat
		r_Earth_Sun = np.linalg.norm(Sun_posECI)
		Rsat = np.linalg.norm(xECI[0:3])
		PV = r_Earth_Sun*SpaceCons.RADIUS_EARTH / (SpaceCons.RADIUS_SUN-SpaceCons.RADIUS_EARTH)
		# angles
		alphaUmb = np.arcsin((SpaceCons.RADIUS_SUN-SpaceCons.RADIUS_EARTH)/r_Earth_Sun) 
		alphaPen = np.arcsin((SpaceCons.RADIUS_SUN+SpaceCons.RADIUS_EARTH)/r_Earth_Sun)
		PV1 = SpaceCons.RADIUS_EARTH / np.arcsin(alphaPen)
		# distance criteria
		cosinEpsilon = np.matmul(Sun_posECI, xECI[0:3]) / r_Earth_Sun / Rsat
		epsilon = np.arccos(cosinEpsilon) # [0,pi]
		satHoriz = Rsat * (-cosinEpsilon)
		satVert = Rsat * np.sin(epsilon)
	
		LengUmb = (PV - satHoriz) * np.tan(alphaUmb)
		LengPen = (PV1 + satHoriz) * np.tan(alphaPen)
	
		# Umbra or Penumbra or light
		# light

		# print(satVert, LengPen)
		# pdb.set_trace() # break ponit
		if satVert >= LengPen:
			pass
		  # tempP_sun = P_sun
		# Umbra
		if satVert <= LengUmb:
			tempP_sun = 0
			# time ratio (0, 1)
			eclipse_ratio = 0
		# penumbra
		if (satVert > LengUmb) and (satVert < LengPen):
			x = LengPen - satVert
			y = x / (LengPen - LengUmb) * SpaceCons.RADIUS_SUN # shadow
			tethaTemp = np.arccos((SpaceCons.RADIUS_SUN-y) / SpaceCons.RADIUS_SUN)
			temp = (tethaTemp - (1-y/SpaceCons.RADIUS_SUN)*np.sin(tethaTemp))/np.pi

			tempP_sun = SpaceCons.P_SUN * (1 - temp)
			# time ratio (0, 1)
			eclipse_ratio = 1 - temp


	### m/s^2
	SRP = 2.0 * tempP_sun * SpaceCons.area2mass * cosinAlpha**2 * vec_n

	return SRP, eclipse_ratio


def SRP_real(xECI, sailControlAngles, Tut1, Ttdb):
	"""
	Calculate SRP due to realistic solar sail in ECI
		inputs:
				xECI(1x6 numpy array) : satellite position & velocity vector in ECI [rx, ry, rz, vx, vy, vz]
				in m
				Tut1, Ttdb(float): 
		outputs:
				SRP_real_n(1x3 numpy array): SRP in sail normal direction, ECI [ai, aj, ak]
				SRP_real_t(1x3 numpy array): SRP in sail transverse direction, ECI [ai, aj, ak]
				SRP: a_SRP_n + a_SRP_t
				eclipse_ratio(float): eclipse condition, 1:light, 0:umbra, (0,1): punumbra
		Ref: Vallado, V4, P574 
			 Mei H, Damaren C J, Zhan X. 
			 End-of-Life Geostationary Satellite Removal Using Realistic Flat Solar Sails.
			 Aerospace Systems, February 2021.
	"""
	# sun position vector in ECI
	Sun_posECI  = thirdBodyPosECI.SunPosECI(Tut1, Ttdb)

	# sunline vector in ECI
	r_sat_Sun = Sun_posECI - xECI[0:3]
	# opposite direction of r_sat_Sun, and normalize
	vec_u = - r_sat_Sun / np.linalg.norm(r_sat_Sun) 


	### sail orientation (control)
	# Construct frame Fs (see ref paper)
	S1 = vec_u
	S3 = (SpaceCons.G3-np.matmul(SpaceCons.G3,vec_u)*vec_u) / np.linalg.norm(SpaceCons.G3-np.matmul(SpaceCons.G3,vec_u)*vec_u)
	S2 = np.cross(S3,S1) / np.linalg.norm(np.cross(S3,S1))
	# Cgs : rotation matrix from Fs to ECI
	Cgs = np.array([[np.matmul(SpaceCons.G1,S1), np.matmul(SpaceCons.G1,S2), np.matmul(SpaceCons.G1,S3)],
					[np.matmul(SpaceCons.G2,S1), np.matmul(SpaceCons.G2,S2), np.matmul(SpaceCons.G2,S3)],
					[np.matmul(SpaceCons.G3,S1), np.matmul(SpaceCons.G3,S2), np.matmul(SpaceCons.G3,S3)]])
	# coordinates of sail normal vector in Fs
	alpha, delta = sailControlAngles[0], sailControlAngles[1]
	ns = np.array([np.cos(alpha),
				   np.sin(alpha) * np.sin(delta),
				   np.sin(alpha) * np.cos(delta)])
	# coordinates of sail normal vector in ECI
	vec_n = np.matmul(Cgs, ns)


	### SRP
	cosinAlpha = np.matmul(vec_n, vec_u)
	# transverse vecter t
	temp_t = vec_u - np.matmul(vec_u, vec_n) * vec_n
	# print(temp_t, np.linalg.norm(temp_t))
	# pdb.set_trace() # break ponit
	vec_t = temp_t / np.linalg.norm(temp_t)
	print('vec_t', vec_t)
	pdb.set_trace() # break ponit

	tempP_sun = SpaceCons.P_SUN
	## eclipse (see ref)
	# eclipse_ratio: 1:light, 0:umbra, (0,1): punumbra
	eclipse_ratio = 1
	if np.matmul(Sun_posECI, xECI[0:3]) < 0: # when the sat is behind the Earth
		# distances: Earth2Sun, Earth2Sat
		r_Earth_Sun = np.linalg.norm(Sun_posECI)
		Rsat = np.linalg.norm(xECI[0:3])
		PV = r_Earth_Sun*SpaceCons.RADIUS_EARTH / (SpaceCons.RADIUS_SUN-SpaceCons.RADIUS_EARTH)
		# angles
		alphaUmb = np.arcsin((SpaceCons.RADIUS_SUN-SpaceCons.RADIUS_EARTH)/r_Earth_Sun) 
		alphaPen = np.arcsin((SpaceCons.RADIUS_SUN+SpaceCons.RADIUS_EARTH)/r_Earth_Sun)
		PV1 = SpaceCons.RADIUS_EARTH / np.arcsin(alphaPen)
		# distance criteria
		cosinEpsilon = np.matmul(Sun_posECI, xECI[0:3]) / r_Earth_Sun / Rsat
		epsilon = np.arccos(cosinEpsilon) # [0,pi]
		satHoriz = Rsat * (-cosinEpsilon)
		satVert = Rsat * np.sin(epsilon)
	
		LengUmb = (PV - satHoriz) * np.tan(alphaUmb)
		LengPen = (PV1 + satHoriz) * np.tan(alphaPen)
	
		# Umbra or Penumbra or light
		# light

		# print(satVert, LengPen)
		# pdb.set_trace() # break ponit
		if satVert >= LengPen:
			pass
		  # tempP_sun = P_sun
		# Umbra
		if satVert <= LengUmb:
			tempP_sun = 0
			# time ratio (0, 1)
			eclipse_ratio = 0
		# penumbra
		if (satVert > LengUmb) and (satVert < LengPen):
			x = LengPen - satVert
			y = x / (LengPen - LengUmb) * SpaceCons.RADIUS_SUN # shadow
			tethaTemp = np.arccos((SpaceCons.RADIUS_SUN-y) / SpaceCons.RADIUS_SUN)
			temp = (tethaTemp - (1-y/SpaceCons.RADIUS_SUN)*np.sin(tethaTemp))/np.pi

			tempP_sun = SpaceCons.P_SUN * (1 - temp)
			# time ratio (0, 1)
			eclipse_ratio = 1 - temp


	### m/s^2
	# SRP in sail normal direction a_n, m^2/s
	temp_a_SRP_n = (1 + SpaceCons.REALSAIL_R_TILD * SpaceCons.REALSAIL_S_TILD) * cosinAlpha**2 \
				   + SpaceCons.REALSAIL_BF * (1 - SpaceCons.REALSAIL_S_TILD) * SpaceCons.REALSAIL_R_TILD * cosinAlpha \
				   + (1 - SpaceCons.REALSAIL_R_TILD) * (SpaceCons.REALSAIL_EPSTONF * SpaceCons.REALSAIL_BF - SpaceCons.REALSAIL_EPSTONB * SpaceCons.REALSAIL_BB) \
				   / (SpaceCons.REALSAIL_EPSTONF + SpaceCons.REALSAIL_EPSTONB) * cosinAlpha
	a_SRP_n = SpaceCons.REFLECTABILITY * tempP_sun * SpaceCons.area2mass * temp_a_SRP_n * vec_n
	 
	# SRP in sail transverse direction a_t, m^2/s
	if cosinAlpha < 1.0:
		# alpha > 0 deg: exists SRP in sail transverse direction
		sinAlpha = (1 - cosinAlpha**2)**(1/2)
		a_SRP_t = SpaceCons.REFLECTABILITY * tempP_sun * SpaceCons.area2mass \
				  * (1 - SpaceCons.REALSAIL_R_TILD * SpaceCons.REALSAIL_S_TILD) \
				  * cosinAlpha * sinAlpha * vec_t
	else:
		a_SRP_t = np.array([0, 0, 0]) # alpha = 0: a_SRP_t = 0
	
	SRP = a_SRP_n + a_SRP_t

	# print('cosinAlpha', cosinAlpha)
	# print('SRP_n, SRP_t:', a_SRP_n, a_SRP_t)
	# pdb.set_trace() # break ponit

	return SRP, eclipse_ratio











