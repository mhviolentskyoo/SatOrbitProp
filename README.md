# SatOrbitProp
Satellite orbit propagator, including:
Control: 
        Control angles (cone and clock angles) of an ideal solar sail 
Disturbances:
        Earth's gravitational disturbances:
				J2 J22 J3 J31 J32 J33
				Thirdbody gravitational accelerations:
				from the Sun and the Moon		

The main file: orbit_propagation.py
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
		satPosECI(1x3 numpy array) : the satellite position vector in ECI [rx,ry,rz]
		in m
		satVelECI(1x3 numpy array) : the satellite velocity vector in ECI [vx,vy,vz]
		in m/s
		startTime(1x6 numpy array) : UTC time [year, mon, day, hr, min, sec]
		timeLength(float)		       : simulation time length
		in sec
		timeStep(float)		         : simulation time step (default = 30)
		in sec
		plotFlag(logical)		       : plot figures or not, True: plot
Returns:
		pos_array(1x3xn numpy array)    : position in ECI [px, py, pz]
		in m
		vel_array(1x3xn numpy array)    : velocity in ECI [vx, vy, vz]
		in m/s
		coe_array(1x3xn numpy array)	  : classical orbital elements [a, ecc, incl, oemga, RAAN, theta]
		epoch_array(floatxn)        	  : simulation time array
		angles_array(1x2xn numpy array) : solar sail control angles
		eclipse_Timeratio(floatxn)      : eclipse time ratio


Reference papers:
[1] Mei H, Damaren C J, Zhan X. Hybrid removal of end-of-life geosynchronous satellites using solar radiation pressure and impulsive thrusts[J]. Advances in Space Research, 2020, 66(4): 974-991.

[2] Mei H, Damaren C J, Zhan X. End-of-life geostationary satellite removal using realistic flat solar sails[J]. Aerospace Systems, 2021: 1-12.
