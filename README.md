# SatOrbitProp

Satellite orbit propagator, used for solar sailing in Earth orbits.

Control: 

    Control angles (cone and clock angles) of an ideal solar sail. (Realistic solar sail parts are not completed.)
        
Disturbances:

    Earth's gravitational disturbances: 
    J2 J22 J3 J31 J32 J33 terms
        
    Third-body gravitational accelerations:
    from the Sun and the Moon		


The main file: orbit_propagation.py

Inputs: 

    Satellite position & velocity in ECI;
  
    Start UTC time [year, mon, day, hr, min, sec];
  
    Simulation time;
    
    Simulation time step.
      
Outputs:

    Satellite orbit propagation using RK4.
      
Notes:

Control: 

    For SRP: 1) Define spacecraft A/m in Spacecons.py. 2) Define solar sail control angles in function: orbit_disturbance.controlAngles.
        
Disturbances:

    1) Earth's gravitational accelerations: J2 J22 J3 J31 J32 J33 terms. 2) Thirdbody gravitational accelerations: from the Sun and the Moon.
    
Reference papers:

[1] Mei H, Damaren C J, Zhan X. Hybrid removal of end-of-life geosynchronous satellites using solar radiation pressure and impulsive thrusts[J]. Advances in Space Research, 2020, 66(4): 974-991.

[2] Mei H, Damaren C J, Zhan X. End-of-life geostationary satellite removal using realistic flat solar sails[J]. Aerospace Systems, 2021: 1-12.

Reference books:

Vallado, V4 & Vallado matlab functions
