#!/usr/bin/env python

import math
import ephem
import sys
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import vec2orbElem
import utils

rad_to_deg = 180.0/math.pi

lsst = ephem.Observer()
lsst.lat = "-39.2333"
lsst.lon = "-70.7333"
lsst.elevation = 2715.0

def testWithMars():
    m = ephem.Mars()

    observations = []
    observations.append(ephem.Date("2014/01/09 3:00"))
    observations.append(observations[-1] + ephem.hour)
    observations.append(ephem.Date("2014/01/11 3:00"))
    observations.append(observations[-1] + ephem.hour)
    observations.append(ephem.Date("2014/01/13 3:00"))
    observations.append(observations[-1] + ephem.hour)

    positions = np.zeros((len(observations), 2))
    absolute_positions = np.zeros((len(observations), 3))
    
    for i in range(len(observations)):
        lsst.date = observations[i]
        m.compute(lsst)
        positions[i, 0] = float(m.ra) * rad_to_deg
        positions[i, 1] = float(m.dec) * rad_to_deg

        print m.hlon, m.hlat, m.ra, m.dec

        # Compute heliocentric position
        absolute_positions[i, 0] = m.sun_distance * np.cos(m.hlon)
        absolute_positions[i, 1] = m.sun_distance * np.sin(m.hlon)
        absolute_positions[i, 2] = m.sun_distance * m.hlat

    print absolute_positions   

    ra_velocities, dec_velocities = getAngularVelocity(observations, positions)
    

    #vec2orbElem(rs, vs, mus)


def getAngularVelocity(observations, positions):

    delta_t = np.zeros(len(observations)/2)
    separations = np.zeros((len(observations)/2, 2))
    angular_velocities = np.zeros(len(observations)/2)

    for i in range(1, len(observations), 2):
        delta_t[(i-1)/2] = (observations[i] - observations[i-1])
        separations[(i-1)/2, 0] = positions[i, 0] - positions[i-1, 0]
        separations[(i-1)/2, 1] = positions[i, 1] - positions[i-1, 1]
    print delta_t
    print separations

    ra_velocities = np.divide(separations[:, 0], delta_t)
    dec_velocities = np.divide(separations[:, 1], delta_t)
    print ra_velocities
    print dec_velocities
    
    fig = plt.figure()
    observed_positions = fig.add_subplot(111)
    observed_positions.plot(positions[:, 0], positions[:, 1], 'bo')
    plt.show()

    return ra_velocities, dec_velocities


def displayFits(observations, positions, elements):
    fig = plt.figure()
    tracks = fig.add_subplot(111)
    tracks.plot(positions[:, 0], positions[:, 1], 'bo')

    time = np.arange(observations[0], observations[-1], 100)

    #for body in elements:
        


def main(argv):
    testWithMars()


if __name__ == '__main__':
    main(sys.argv[1:])
