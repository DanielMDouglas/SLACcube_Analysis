clock_interval = 0.1 # us/tick -- 10 MHz clock rate
v_drift = 1.648 # mm/us (rough estimate)

detector_bounds = [[-150, 150], [-150, 150], [0, 300]] # mm (x, y, z)

drift_distance = detector_bounds[2][1] - detector_bounds[2][0] 

drift_window = drift_distance/(v_drift*clock_interval) # maximum drift time

drift_direction = 1 # +/- 1 depending on the direction of the drift in z
