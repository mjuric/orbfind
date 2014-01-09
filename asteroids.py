#!/usr/bin/env python

import math
import ephem
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

rad_to_deg = 180.0/math.pi

lsst = ephem.Observer()
lsst.lat = "-39.2333"
lsst.lon = "-70.7333"
lsst.elevation = 2715.0

observations = []
observations.append(ephem.Date("2014/01/09 3:00"))
observations.append(observations[-1] + ephem.hour)
observations.append(ephem.Date("2014/01/11 3:00"))
observations.append(observations[-1] + ephem.hour)
observations.append(ephem.Date("2014/01/13 3:00"))
observations.append(observations[-1] + ephem.hour)

delta_t = observations[1:5] - observations[0:4]
print delta_t

for date in observations:
    print 

m = ephem.Mars()
positions = []

for obs in observations:
    lsst.date = obs
    m.compute(lsst)
    positions.append((float(m.ra) * rad_to_deg, float(m.dec)*rad_to_deg))

ra = np.array([pos[0] for pos in positions])
dec = np.array([pos[1] for pos in positions])
print positions
fig = plt.figure()
observed_positions = fig.add_subplot(111)
observed_positions.plot(ra, dec, 'bo')
plt.show()
