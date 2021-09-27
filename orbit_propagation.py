import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pdb

import orbit_disturbance
import SpaceCons
import SpaceConvert


class satOrbit_PROP:
	"""
	A class defined for orbit propagation.
	Input : 
			satellite position & velocity in ECI
			start UTC time [year, mon, day, hr, min, sec]
			simulation time
			simulation time step
	Output:
			propagate satellite orbit using RK4
	Note:
		Control: 
				SRP
				Define spacecraft A/m in Spacecons.py
				Define solar sail control angles in function: orbit_disturbance.controlAngles
		Disturbances including:
				Earth's gravitational accelerations:
				J2 J22 J3 J31 J32 J33
				Thirdbody gravitational accelerations:
				from the Sun and the Moon			
	Args:
		satPosECI(1x3 numpy array): the satellite position vector in ECI [rx,ry,rz]
		in m
		satVelECI(1x3 numpy array): the satellite velocity vector in ECI [vx,vy,vz]
		in m/s
		startTime(1x6 numpy array): UTC time [year, mon, day, hr, min, sec]
		timeLength(float)		  : simulation time length
		in sec
		timeStep(float)		      : simulation time step (default = 30)
		in sec
		plotFlag(logical)		  : plot figures or not, True: plot
	Returns:
		pos_array(1x3xn numpy array)    : position in ECI [px, py, pz]
		in m
		vel_array(1x3xn numpy array)    : velocity in ECI [vx, vy, vz]
		in m/s
		coe_array(1x3xn numpy array)	: classical orbital elements [a, ecc, incl, oemga, RAAN, theta]
		epoch_array(floatxn)        	: simulation time array
		angles_array(1x2xn numpy array) : solar sail control angles
		eclipse_Timeratio(floatxn)      : eclipse time ratio
	"""

	def __init__(self, satPosECI, satVelECI, startTime, 
				 timeLength, timeStep=30, plotFlag='True'):
		"""
		Initiate
		Inputs:
				satPosECI(1x3 numpy array)   : satellite initial position in ECI [px, py, pz]
				satVelECI(1x3 numpy array)   : satellite initial velocity in ECI [vx, vy, vz]
				startTime(1x6 numpy array)   : start UTC time
				timeLength(float)            : simulation time length
				timeStep(float)  			 : time step in numerical integrator
				plotFlag(boolean)			 : plot figs or not
		"""
		self._position_initial = satPosECI
		self._velocity_initial = satVelECI
		self._t_start  		   = startTime
		self._t_sim 		   = timeLength
		self._h 			   = timeStep
		self._plotFlag		   = plotFlag
 
 	# orbit propagate
	def orbit_prop(self):
		"""
		Numerical orbit propagation using RK4
		"""
		# initiate
		print('Simulation time length:', self._t_sim/3600/24, 'days')
		xECI_initial = np.array([*self._position_initial, *self._velocity_initial])

		# numeriacl integration
		pos_array, vel_array, coe_array, epoch_array, angles_array, eclipse_Timeratio \
			= self.rk4(xECI_initial, self._t_start, self._t_sim, self._h)

		return pos_array, vel_array, coe_array, epoch_array, angles_array, eclipse_Timeratio

	# plot figs
	def plot_figs(self):
		if self._plotFlag == 'True':
			### plot orbital elements
			# semi-major axis
			plt.figure() 
			plt.plot(epoch_array / 3600 / 24, (coe_array[:,0] - SpaceCons.RGEO) / 1e3,
					 linestyle='-.', label='semi-major axis vs time') 
					 # time: number of days
			plt.legend(loc='upper right')
			plt.xlabel('time (days)')
			plt.ylabel('semi-major axis - RGEO (km)')
			# plt.xlim((-0.5, 21.5))
			# plt.ylim((-1, 3))
			plt.show()

			# eccentricity
			plt.figure() 
			plt.plot(epoch_array / 3600 / 24, coe_array[:,1],
					 linestyle='-.', label='eccentricity vs time')
			plt.legend(loc='upper right')
			plt.xlabel('time (days)')
			plt.ylabel('eccentricity')
			# plt.xlim((-0.5, 21.5))
			# plt.ylim((-0.1, 0.1))
			plt.show()

			### plot 3D trajectory
			fig = plt.figure()
			ax  = plt.axes(projection ='3d')

			ax.plot3D(pos_array[:,0], pos_array[:,1], pos_array[:,2],
					  linestyle ='-.', label ='3D trajectory')
			ax.set_xlabel('xECI')
			ax.set_ylabel('yECI')
			ax.set_zlabel('zECI')
			plt.show()

			### print terminal states
			# print terminal ECI coordinates
			print("\nTerminal position in ECI (in m):", 
					pos_array[pos_array.shape[0]-1, 0],
					pos_array[pos_array.shape[0]-1, 1],
					pos_array[pos_array.shape[0]-1, 2])
			# print("\nTerminal orbital elements:", 
			# 	  "\nSemimajor axis (m):",  	  coe_array[coe_array.shape[0]-1, 0],
			# 	  "\nEccentricity :",	   	  coe_array[coe_array.shape[0]-1, 1],
			# 	  "\nInclination (deg):",		 coe_array[coe_array.shape[0]-1, 2]*180/np.pi,
			# 	  "\nArgument of parigee (deg):", coe_array[coe_array.shape[0]-1, 3]*180/np.pi,
			# 	  "\nRAAN (deg):",				coe_array[coe_array.shape[0]-1, 4]*180/np.pi,
			# 	  "\nTrue anomaly (deg):",		coe_array[coe_array.shape[0]-1, 5]*180/np.pi)
			# terminal orbital elements
			print("\nTerminal orbital elements:", 
				  "\nSemimajor axis (m):",  	  coe_array[coe_array.shape[0]-1, 0],
				  "\nEccentricity :",	   	  coe_array[coe_array.shape[0]-1, 1],
				  "\nInclination (rad):",		 coe_array[coe_array.shape[0]-1, 2],
				  "\nArgument of parigee (rad):", coe_array[coe_array.shape[0]-1, 3],
				  "\nRAAN (rad):",				coe_array[coe_array.shape[0]-1, 4],
				  "\nTrue anomaly (rad):",		coe_array[coe_array.shape[0]-1, 5])

			### plot control angles
			# alpha, deg
			plt.figure() 
			plt.plot(epoch_array / 3600 / 24, angles_array[:,0]*180/np.pi,
					 linestyle='-.', label=r'$ \alpha $')
			plt.legend(loc='upper left')
			plt.xlabel('time (days)')
			plt.ylabel(r'$ \alpha $')
			# plt.xlim((-0.5, 21.5))
			# plt.ylim((-0.1, 0.1))
			plt.show()
			# beta, deg
			plt.figure() 
			plt.plot(epoch_array / 3600 / 24, angles_array[:,1]*180/np.pi,
					 linestyle='-.', label=r'$ \beta $')
			plt.legend(loc='upper left')
			plt.xlabel('time (days)')
			plt.ylabel(r'$ \beta $')
			# plt.xlim((-0.5, 21.5))
			# plt.ylim((-0.1, 0.1))
			plt.show()



	""" Numerical Integrators """
	# RK4, Runge-Kutta 4th Order Numerical Integrator
	def rk4(self, s, t_start, t_sim, h=30):
		"""
		Runge-Kutta 4th Order Numerical Integrator
		Args:
			s(1x6 numpy array): the initial state vector [rx,ry,rz,vx,vy,vz]
			t_start(float) : (1x6 numpy array): [year, mon, day, hr, min, sec]
			t_sim(float)   : simulation time, sec
			h(float)	   : step-size
		Returns:
			1x6 numpy array: the state at time tf
		"""

		pos, vel, coe, epoch, controlAngles, eclipse_Timeratio = [],[],[],[],[],[]

		t0 = t_start[5] # sec
		tf = t0 + t_sim # final sec

		t_sec = t0

		# orbit propagation
		while(abs(tf-t_sec) > 0.00001):
			# if tf - t < h, propagate time reminder
			if (abs(tf-t_sec) < abs(h)):
				h = tf-t_sec

			# actual time
			t_current = np.array([t_start[0],
								  t_start[1],
								  t_start[2],
								  t_start[3],
								  t_start[4],
								  t_sec])
			t_half_h  = np.array([t_start[0],
								  t_start[1],
								  t_start[2],
								  t_start[3],
								  t_start[4],
								  t_sec + h/2])
			t_h		  = np.array([t_start[0],
								  t_start[1],
								  t_start[2],
								  t_start[3],
								  t_start[4],
								  t_sec + h])

			k1, SRP_type = self.state_derivative(s,		     t_current)
			k2, __	     = self.state_derivative(s + h*k1/2, t_half_h)
			k3, __ 		 = self.state_derivative(s + h*k2/2, t_half_h)
			k4, __   	 = self.state_derivative(s + h*k3,   t_h)

			# save position & velocity vectors
			# before state update
			pos.append(s[0:3])
			vel.append(s[3:6])

			# convert to classical orbital elements and save
			# before state update
			coe_full = SpaceConvert.RV2OrbitalEle(s)
			coe6 = [coe_full.a, 
					coe_full.ecc, 
					coe_full.incl, 
					coe_full.omega, 
					coe_full.RAAN, 
					coe_full.theta]
			coe.append(coe6)
			# print("coe6", coe6)
			# pdb.set_trace() # break ponit

			# save time epochs & eclipse time ratio
			# before update
			epoch.append(t_sec)

			controlAngles.append(SRP_type[0:2])
			eclipse_Timeratio.append(SRP_type[2])
			# print("SRP_type[0:2]", SRP_type[0:2])
			# pdb.set_trace() # break ponit

			### orbit propagation
			s	  = s + h/6 * (k1+2*k2+2*k3+k4)
			t_sec = t_sec+h

		### save pos,vel,coe,epoch to array
		pos_array   	  = np.asarray(pos)
		vel_array   	  = np.asarray(vel)
		coe_array	 	  = np.asarray(coe)
		epoch_array   	  = np.asarray(epoch)
		angles_array	  = np.asarray(controlAngles)
		eclipse_Timeratio = np.asarray(eclipse_Timeratio)

		return pos_array, vel_array, coe_array, epoch_array, angles_array, eclipse_Timeratio

	# RKF45, Runge-Kutta Fehlberg 4(5) Numerical Integrator
	# not completed
	def rkf45(self, s, t0, tf, h=10, tol=1e-6):
		"""
		********** Not completed **********
		Runge-Kutta Fehlberg 4(5) Numerical Integrator
		Args:
			s(1x6 numpy array): the initial state vector [rx,ry,rz,vx,vy,vz]
			t0(float)  : initial time
			tf(float)  : final time
			h(float)   : step-size
			tol(float) : tolerance of error
		Returns:
			1x6 numpy array: the state at time tf
		"""

		t = t0

		# backward propagation
		if tf < t0: 
			h = -h

		# orbit propagation
		while(abs(tf-t) > 0.00001):
			# if tf - t < h, propagate time reminder
			if (abs(tf-t) < abs(h)):
				h = tf-t

			k1 = h*state_derivative(s)
			k2 = h*state_derivative(s+k1/4)
			k3 = h*state_derivative(s+3/32*k1+9/32*k2)
			k4 = h*state_derivative(s+1932/2197*k1-7200/2197*k2+7296/2197*k3)
			k5 = h*state_derivative(s+439/216*k1-8*k2+3680/513*k3-845/4104*k4)
			k6 = h*state_derivative(s-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)

			y = s+25/216*k1+1408/2565*k3+2197/4104*k4-k5/5
			z = s+16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6

			s = y
			t = t+h

			err = np.linalg.norm(y-z)
			h = h*0.84*(tol/err)**0.25

		return s



	""" Find the state derivative """
	def state_derivative(self, xECI, time):
		"""
		Returns the time derivative of a given state 
		Args:
			xECI(1x6 numpy array): the state vector in ECI [rx,ry,rz,vx,vy,vz]
			in m, m/s
			time(1x6 numpy array): (1x6 numpy array): [year, mon, day, hr, min, sec]
		Returns:
			1x6 numpy array: the time derivative of s [vx,vy,vz,ax,ay,az]
			in m/s, m/s^2
		"""

		# orbit radius, in m 
		r = np.linalg.norm(xECI[0:3])

		### time convert
		# find Tut1, Ttdb
		Tut1, Ttdb = SpaceConvert.TimeConvert.convtime(time[0],  # year
													   time[1],  # mon
													   time[2],  # day
													   time[3],  # hr
													   time[4],  # min
													   time[5])  # sec
		# find GMST, in rad
		GMST = SpaceConvert.TimeConvert.Tut1_2_GMST(Tut1)


		### Earth's gravitational acceleration
		# find pos in ECEF
		posECEF = SpaceConvert.ECI2ECEF(xECI[0:3], GMST)

		# twobody acceleration in ECEF, in m/s
		a_twoBody = -SpaceCons.MU_EARTH / (r**3) * posECEF

		# perturbative accelerations in ECEF
		a_J2  = orbit_disturbance.J2_pert(posECEF)
		a_J22 = orbit_disturbance.J22_pert(posECEF)
		a_J3  = orbit_disturbance.J3_pert(posECEF)
		a_J31 = orbit_disturbance.J31_pert(posECEF)
		a_J32 = orbit_disturbance.J32_pert(posECEF)
		a_J33 = orbit_disturbance.J33_pert(posECEF)

		# total gravitational acceleration in ECEF
		# a = a_twoBody # tow-body problem, debug: correct
		aECEF = a_twoBody + a_J2 + a_J22 + a_J3 + a_J31 + a_J32 + a_J33
		# convert accelerations to ECI
		aEarthECI = SpaceConvert.ECEF2ECI(aECEF, GMST)

		# print("a_ECEF", aECEF)
		# print("aECI", aECI)
		# pdb.set_trace() # break ponit


		### Thirdbody's gravitational acceleration
		# Sun & Moon position vectors in ECI
		Sun_posECI  = orbit_disturbance.thirdBodyPosECI.SunPosECI(Tut1, Ttdb)
		Moon_posECI = orbit_disturbance.thirdBodyPosECI.MoonPosECI(Ttdb)
		# Thirdbody gravitational acceleration vectors in ECI
		acc_GraSun  = orbit_disturbance.thirdBodyAccECI.thirdbody_Sun(xECI[0:3], 
																	  Sun_posECI)
		acc_GraMoon = orbit_disturbance.thirdBodyAccECI.thirdbody_Moon(xECI[0:3], 
																	   Moon_posECI)
		# Thirdbody's gravitational acceleration
		aThirdbodyECI = acc_GraSun + acc_GraMoon

		# print('Sun_posECI', Sun_posECI)
		# print('Moon_posECI', Moon_posECI)
		# print('acc_GraSun', acc_GraSun)
		# print('acc_GraMoon', acc_GraMoon)
		# pdb.set_trace() # break ponit


		### Solar radiation pressure (from a solar sail)
		# solar sail control angles
		sailControlAngles = orbit_disturbance.controlAngles(xECI, Tut1, Ttdb)
		# SRP for an ideal solar sail
		aSRPECI, eclipse_ratio = orbit_disturbance.SRP_perfect(xECI, sailControlAngles, Tut1, Ttdb)
		# print('aSRPECI, eclipse_ratio', aSRPECI, eclipse_ratio)
		# pdb.set_trace() # break ponit


		### state derivative
		# total acceleration in ECI
		aECI = aEarthECI + aThirdbodyECI + aSRPECI
		# save states
		xECI_dot = np.array([*xECI[3:6], *aECI]) # *把两个[][]去掉了
		# save sailControlAngles and eclipse_ratio
		SRP_type = np.array([*sailControlAngles, eclipse_ratio])

		return xECI_dot, SRP_type



""" main """
if __name__ == "__main__":

	# initial states, in m
	# e.g. a GEO satellite
	GEO_position_initial = np.array([ 0.0, 
									  SpaceCons.RGEO, 
									  1])
	GEO_velocity_initial = np.array([-SpaceCons.OMEGA_EARTH*SpaceCons.RGEO, 
									 -SpaceCons.OMEGA_EARTH*0.0, 
									  0])
	# # initial state
	# xECI_initial = np.array([*position_initial, *velocity_initial])

	# initial orbital elements
	# coe_ini = SpaceConvert.RV2OrbitalEle(xECI_initial)
	# print(coe_ini.a, coe_ini.ecc, coe_ini.incl, coe_ini.omega, coe_ini.RAAN, coe_ini.theta)

	# initial time, (1x6 numpy array): [year, mon, day, hr, min, sec]
	t_start = np.array([2017, 1, 1, 0, 0, 0])
	# simulation time, sec
	t_sim = 3600 * 24 * 1 # 3600 * 24 * days
	# print('Simulation time length:', t_sim/3600/24, 'days')
	# time step
	h  = 30 # sec
	plot_flag = 'True'

	GEOSat_prop = satOrbit_PROP(GEO_position_initial,
								GEO_velocity_initial,
								t_start,
								t_sim,
								h,
								plot_flag)

	pos_array, vel_array, coe_array, epoch_array, angles_array, eclipse_Timeratio \
			= GEOSat_prop.orbit_prop()

	GEOSat_prop.plot_figs()

	






























