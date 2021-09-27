"""
Time & Frame convert functions, includes:
RV2OrbitalEle(s): 
				   R & V in ECI to classical orbital elements
Rotation Matrices: 
				   class RotMat
				   def Cx(angle) Cy Cz
ECI2ECEF(XECI, GMST) : get XECI coordinates in ECEF
"""

import numpy as np
import math
import SpaceCons


""" R & V in ECI to classical orbital elements """
def RV2OrbitalEle(s):
	"""
	inputs: 
			s(1x6 numpy array): the state vector [rx,ry,rz,vx,vy,vz]
			in m
	outputs: 
			structure COE containing:
			a , ecc, incl, epsilon, RAAN, omega, nj, ek, theta
	"""

	class COE:
		def __init__(self, a , ecc, incl, epsilon, RAAN, omega, nj, ek, theta):
			self.a	   = a
			self.ecc	 = ecc
			self.incl	= incl
			self.epsilon = epsilon
			self.RAAN	= RAAN
			self.omega   = omega
			self.nj	  = nj
			self.ek	  = ek
			self.theta   = theta

	# orbit radius, in m
	R = np.linalg.norm(s[0:3])

	# orbital energy
	epsilon = np.matmul(s[3:6],s[3:6]) / 2 - SpaceCons.MU_EARTH / R
	COE.epsilon = epsilon
	
	# semi-major axis, in m
	COE.a = - SpaceCons.MU_EARTH / (2 * epsilon)

	# eccentricity
	hHat = np.cross(s[0:3], s[3:6])
	eHat = np.cross(s[3:6],hHat) / SpaceCons.MU_EARTH - s[0:3]/R
	ecc  = np.linalg.norm(eHat)
	COE.ecc = ecc
	
	# inclination, in rad
	h = np.linalg.norm(hHat)
	COE.incl = np.arccos(hHat[2] / h)

	# RAAN, in rad
	nHat = np.cross(SpaceCons.G3, hHat)
	n = np.linalg.norm(nHat)
	cosRAAN = nHat[0] / n
	RAAN = np.arccos(cosRAAN)
	if  nHat[1] < 0 :
		RAAN = 2 * math.pi - RAAN
	COE.RAAN = RAAN

	# omega (arguement of parigee): in rad
	cosOmega = np.matmul(nHat, eHat) / n / ecc
	omega = np.arccos(cosOmega)
	if eHat[2] < 0 :
		omega = 2 * math.pi - omega
	COE.omega = omega

	# ture anomaly, in rad
	cosTheta = np.matmul(s[0:3], eHat) / ecc / R
	theta = np.arccos(cosTheta)
	if np.matmul(s[0:3], s[3:6]) < 0 :
		theta = 2 * math.pi - theta
	COE.theta = theta

	# nj & ek
	COE.nj = nHat[1]
	COE.ek = eHat[2]

	return COE



""" Rotation Matrices """
class RotMat:

	def __init__(self):
		pass

	@staticmethod
	def Cx(angle):
		# angle, in rad
		RotMatCx = np.array([[1,  0,			  0],
							 [0,  np.cos(angle),  np.sin(angle)],
					  		 [0, -np.sin(angle),  np.cos(angle)]])
		return RotMatCx

	@staticmethod
	def Cy(angle):
		# angle, in rad
		RotMatCy = np.array([[np.cos(angle),  0,  -np.sin(angle)],
							 [0,			  1,   0],
							 [np.sin(angle),  0,   np.cos(angle)]])
		return RotMatCy

	@staticmethod
	def Cz(angle):
		# angle, in rad
		RotMatCz = np.array([[ np.cos(angle),  np.sin(angle),  0],
							 [-np.sin(angle),  np.cos(angle),  0],
							 [ 0,			   0,			   1]])
		return RotMatCz


""" ECI2ECEF """
def ECI2ECEF(XECI, GMST):
	"""
	inputs: 
			XECI(1x3 numpy array): coordinates in ECI [x,y,z]
			in m
			GMST: Greenwich mean siderial time
			in rad
	outputs: 
			XECEF(1x3 numpy array): coordinates in ECEF [x,y,z]
			in m
	Notes: 
			It's a simplified model
	"""

	# rotation matrix Cz, GMST in rad
	rot_mat = RotMat.Cz(GMST)

	XECEF = np.matmul(rot_mat, XECI)

	return XECEF

""" ECEF2ECI """
def ECEF2ECI(XECEF, GMST):
	"""
	inputs: 
			XECI(1x3 numpy array): coordinates in ECI [x,y,z]
			in m
			GMST: Greenwich mean siderial time
			in rad
	outputs: 
			XECEF(1x3 numpy array): coordinates in ECEF [x,y,z]
			in m
	Notes: 
			It's a simplified model
	"""

	# rotation matrix Cz, GMST in rad
	rot_mat = RotMat.Cz(-GMST)

	XECI = np.matmul(rot_mat, XECEF)

	return XECI



""" Time convert functions """
class TimeConvert:

	def __init__(self):
		pass

	@staticmethod
	def hms2sec(hr, minutes, sec):
		"""
			Converts hours, minutes and seconds into seconds 
			from the beginning of the day.
			inputs:
					hr	  - hours	-  0 .. 24
					min   - minutes  -  0 .. 59
					sec   - seconds  -  0.0 .. 59.99
			outputs:
					utsec - seconds  -  0.0 .. 86400.0
		"""
		utsec  = hr * 3600.0 + minutes * 60.0 + sec

		return utsec

	@staticmethod
	def sec2hms(utsec):
		"""
			Converts seconds from the beginning of the day into hours,
			minutes and seconds.
			iputs:
					utsec - seconds  -  0.0 .. 86400.0
			outputs:
					hr	- hours	-  0 .. 24
					min   - minutes  -  0 .. 59
					sec   - seconds  -  0.0 .. 59.99
		"""
		temp    = utsec / 3600.0
		hr      = np.fix( temp )
		minutes = np.fix( (temp - hr)*60.0 )
		sec     = (temp - hr - minutes/60.0) * 3600.0

		return hr, minutes, sec

	@staticmethod
	def jday(yr, mon, day, hr, minutes, sec):
		"""
			Find the julian date given the year, month, day, and time.
			inputs:
					year	- year				   - 1900 .. 2100
					mon	    - month				   - 1 .. 12
					day	    - day				   - 1 .. 28,29,30,31
					hr	    - universal time hour  - 0 .. 23
					min	    - universal time min   - 0 .. 59
					sec	    - universal time sec   - 0.0 .. 59.999
			outputs:
					jd	  - julian date					- days from 4713 bc
					jdfrac  - julian date fraction of a day  - 0.0 to 1.0
			Ref:
					vallado	   2007, 189, alg 14, ex 3-14
					Aerodynamics software by David Vallado
		"""
		jd = 367.0 * yr \
			 - np.floor( (7 * (yr + np.floor( (mon + 9) / 12.0) ) ) * 0.25 ) \
			 + np.floor( 275 * mon / 9.0 ) \
			 + day + 1721013.5 
		jdfrac = (sec + minutes * 60.0 + hr *3600.0) / 86400.0

		# check jdfrac
		if jdfrac > 1.0 :
			jd = jd + np.floor(jdfrac)
			jdfrac = jdfrac - np.floor(jdfrac)

		return jd, jdfrac

	@staticmethod
	def convtime(year, mon, day, hr, minutes, sec,
				 timezone=SpaceCons.TIMEZONE, dut1=SpaceCons.DUT1, dat=SpaceCons.DAT):
		"""
			Finds the time parameters and julian century values for inputs
			of utc time
			inputs: utc time
					year		- year						   1900 .. 2100
					mon		 - month						  1 .. 12
					day		 - day							1 .. 28,29,30,31
					hr		  - universal time hour			0 .. 23
					min		 - universal time min			 0 .. 59
					sec		 - universal time sec (utc)	   0.0  .. 59.999
					timezone	- offset to utc from local site  0 .. 23 hr
					dut1		- delta of ut1 - utc			 sec
					dat		 - delta of tai - utc			 sec
			outputs:
					ut1		 - universal time				 sec
					tut1		- julian centuries of ut1
					jdut1	   - julian date (days only)		days from 4713 bc
					jdut1Frac   - julian date (fraction of a day)days from 0 hr of the day
					utc		 - coordinated universal time	 sec
					tai		 - atomic time					sec
					tdt		 - terrestrial dynamical time	 sec
					ttdt		- julian centuries of tdt
					jdtt		- julian date (days only)		days from 4713 bc
					jdttFrac	- julian date (fraction of a day)days from 0 hr of the day
					tdb		 - terrestrial barycentric time   sec
					ttdb		- julian centuries of tdb
					jdtdb	   - julian date of tdb			 days from 4713 bc
					tcb		 - celestial barycentric time	 sec
					tcg		 - celestial geocentric time	  sec
					jdtdb	   - julian date (days only)		days from 4713 bc
					jdtdbFrac   - julian date (fraction of a day)days from 0 hr of the day
			Ref:
					vallado	   2007, 201, alg 16, ex 3-7
					Aerodynamics software by David Vallado
		"""
		localhr = timezone + hr
		utc = TimeConvert.hms2sec(localhr, minutes, sec)

		ut1 = utc + dut1
		hrtemp,mintemp,sectemp = TimeConvert.sec2hms(ut1)
		JDut1, JDut1frac = TimeConvert.jday(year, mon, day, hrtemp, mintemp, sectemp)
		
		Tut1 = (JDut1 + JDut1frac - 2451545.0) / 36525.0

		tai = utc + dat
		hrtemp,mintemp,sectemp = TimeConvert.sec2hms(tai)
		jdtai, jdtaifrac = TimeConvert.jday(year, mon, day, hrtemp, mintemp, sectemp)

		tt = tai + 32.184
		hrtemp,mintemp,sectemp = TimeConvert.sec2hms(tt)
		JDtt, JDttfrac = TimeConvert.jday(year, mon, day, hrtemp, mintemp, sectemp)
		Ttt = (JDtt + JDttfrac - 2451545.0)/ 36525.0

		tdb = tt + 0.001658 * np.sin(628.3076*Ttt+6.2401)  \
				 + 0.000022 * np.sin(575.3385*Ttt+4.2970)  \
				 + 0.000014 * np.sin(1256.6152*Ttt+6.1969) \
				 + 0.000005 * np.sin(606.9777*Ttt+4.0212)  \
				 + 0.000005 * np.sin(52.9691*Ttt+0.4444)   \
				 + 0.000002 * np.sin(21.3299*Ttt+5.5431)   \
				 + 0.000010 * Ttt * np.sin(628.3076*Ttt+4.2490)
		hrtemp,mintemp,sectemp = TimeConvert.sec2hms(tdb)
		JDtdb, JDtdbfrac = TimeConvert.jday(year, mon, day, hrtemp, mintemp, sectemp)
		Ttdb = (JDtdb + JDtdbfrac - 2451545.0)/ 36525.0

		tcg = tt + 6.969290134e-10*(jdtai - 2443144.5)*86400.0
		hrtemp,mintemp,sectemp = TimeConvert.sec2hms(tcg)
		jdtcg, jdtcgfrac = TimeConvert.jday(year, mon, day, hrtemp, mintemp, sectemp)
		tt2 = tcg - 6.969290134e-10*(jdtcg+jdtcgfrac-2443144.5003725)*86400.0

		tcbmtdb = 1.55051976772e-8*(jdtai+jdtaifrac - 2443144.5)*86400.0
		tcb = tdb + tcbmtdb
		hrtemp,mintemp,sectemp = TimeConvert.sec2hms(tcb)
		jdtcb, jdtcbfrac = TimeConvert.jday(year, mon, day, hrtemp, mintemp, sectemp)

		return Tut1, Ttdb

	@staticmethod
	def Tut1_2_GMST(Tut1):
		"""
			From Tut1 to GMST, GMST in rad
			inputs:
					Tut1
			outputs:
					GMST in rad
			Ref: 
					Vallado, V4, P189
					Vallado matlab functions, gstime & ex3_5
		"""
		# GMST in sec
		GMST_sec = - 6.2e-6 * Tut1 * Tut1 * Tut1 + 0.093104 * Tut1 * Tut1 \
				   + (876600.0 * 3600.0 + 8640184.812866) * Tut1 + 67310.54841
		 
		# GMST in degrees
		temp = np.remainder(GMST_sec, 86400)

		GMST_degree = temp / 240.0 
		 
		# check quadrants
		if GMST_degree < 0.0 :
			 GMST_degree = GMST_degree + 360.0

		GMST = GMST_degree * math.pi / 180

		return GMST



""" for test """
if __name__ == "__main__":

	XECI  = np.array([1,1,3])
	GMST  = math.pi / 2

	XECEF = ECI2ECEF(XECI, GMST)

	print(XECEF)




