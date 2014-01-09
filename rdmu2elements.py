import numpy
import vec2orbElem
import pdb

# to run me: 
# reload(vec2orbElem) ; reload(rdmu2elements) ; 

# this is a first step toward getting all of the 6D orbital elements allowed by the observation of a particular observation of a solar system moving object and its motion

# here we pretend we are observing from the sun so that the mapping of r, d, mu_r, mu_d into position, velocity space is straightforward

# the result is a lot of orbital elements computed over a grid of possible distances and radial velocities

# at this stage you can do
# res = rdmu2elements.rdmu2elements(30, 0, 30, 30)
# to get all of the possible orbital elements corresponding to an observation of an object at r, d = (30, 0), and mu_r, mu_d = (30,30)  I screwed up the units there and the units of mu are senseless and need to be fixed
# orbital elements are from vec2orbElem and ~described there

def rdmu2elements(r, d, mur, mud, ddist=.1, dv=.1):
    auinm = 149597870700.
    
    distrange = [0.01, 20] # AU (~1000 points)
    velrange = [-80,80] # vel (~6000 points)
    ndist = (distrange[1]-distrange[0])/ddist
    nvel = (velrange[1]-velrange[0])/dv
    dist = numpy.linspace(distrange[0], distrange[1], ndist)
#    dist = numpy.array([1])
#    ndist = 1
    vel = numpy.linspace(velrange[0], velrange[1], nvel)
    posg = numpy.zeros((3, ndist, nvel))
    velg = numpy.zeros((3, ndist, nvel))
    rr = numpy.radians(r)
    dr = numpy.radians(d)
    posg[0,:,:] = (numpy.cos(rr)*numpy.cos(dr)*dist).reshape(-1,1)
    posg[1,:,:] = (numpy.sin(rr)*numpy.cos(dr)*dist).reshape(-1,1)
    posg[2,:,:] = (numpy.sin(dr)*dist).reshape(1,-1,1)

    velg[0,:,:] = ((-numpy.sin(rr)*dist*mur*numpy.pi/180.).reshape(-1, 1) +
                   (-numpy.sin(rr)*dist*mud*numpy.pi/180.).reshape(-1, 1) +
                   numpy.cos(rr)*numpy.sin(dr)*vel.reshape(1,-1))
    velg[1,:,:] = ((numpy.cos(rr)*dist*mur*numpy.pi/180.).reshape(-1, 1) +
                   (numpy.sin(rr)*numpy.sin(dr)*dist*mud*numpy.pi/180.).reshape(-1, 1) +
                   numpy.sin(rr)*numpy.cos(dr)*vel.reshape(1,-1))
    velg[2,:,:] = ((numpy.pi/180.*dist*numpy.cos(dr)*mud).reshape(-1, 1) +
                   numpy.sin(dr)*vel.reshape(1,-1))
    return posg, velg, vec2orbElem.vec2orbElem(posg.reshape(3,-1)*auinm, velg.reshape(3,-1)*1000, mus=1.3275e20)
