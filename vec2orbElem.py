import numpy
from numpy import sin, cos
import pdb

def vec2orbElem(rs, vs, mus):
# %vec2orbElem - Convert position and velocity vectors to orbital elements.
# %      vec2orbElem(rs,vs,mus) converts positions (rs) and velocities (vs)
# %      of bodies with gravitational parameters (mus) to Keplerian orbital
# %      elements.
# %
# %      INPUTS:
# %      rs       3n x 1 stacked initial position vectors:
# %               [r1(1);r1(2);r1(3);r2(1);r2(2)r2(3);...;rn(1);rn(2);rn(3)]
# %               or 3 x n matrix of position vecotrs.
# %      vs       3n x 1 stacked initial velocity vectors or 3 x n matrix
# %      mus      gravitational parameters (G*m_i) where G is the
# %               gravitational constant and m_i is the mass of the ith body.
# %               if all vectors represent the same body, mus may be a
# %               scalar.
# %      OUTPUTS:
# %      a        semi-major axes
# %      e        eccentricities
# %      E        eccentric anomalies
# %      I        inclinations
# %      omega    arguments of pericenter
# %      Omega    longitudes of ascending nodes
# %      P        orbital periods
# %      tau      time of periapsis crossing
# %      A, B     orientation matrices (see Vinti, 1998)
# %
# %      All units must be complementary, i.e., if positions are in AUs, and
# %      time is in days, dx0 must be in AU/day, mus must be in
# %      AU^3/day^2 (these are the units in solarSystemData.mat).
# %
# %      The data in solarSystemData.mat was downloaded from JPL's System Web
# %      Interface (http://ssd.jpl.nasa.gov/?horizons).  It includes
# %      positions for the planets, the sun and pluto (because I went to
# %      grade school before 2006).  Positions for planets with moons are for
# %      the barycenters. 
# %
# %      Example:
# %      %solar system oribtal elements
# %      ssdat = load('solarSystemData.mat');
# %      rs = ssdat.p0(1:end-3) - repmat(ssdat.p0(end-2:end),9,1);
# %      vs = ssdat.v0(1:end-3) - repmat(ssdat.v0(end-2:end),9,1);
# %      mus = ssdat.mus(1:9) + ssdat.mus(10);
# %      [a,e,E,I,omega,Omega,P,tau,A,B] = vec2orbElem(rs,vs,mus);
# %      %convert back:
# %      r = A*diag(cos(E) - e) + B*diag(sin(E));
# %      rdot = (-A*diag(sin(E))+B*diag(cos(E)))*...
# %              diag(sqrt(mus(:).'.*a.^-3)./(1 - e.*cos(E)));
# %
# % See also: atan2

# % Written by Dmitry Savransky, 9 July 2008, dsavrans@princeton.edu
# % Update and rewrite w/ switch to atan2, 18 Feb 2011 ds
# hackily ported to python 2014 Jan 9 EFS

# %condition inputs
    nplanets = rs.shape[1]

    v2s = numpy.sum(vs**2., axis=0)
    r = numpy.sqrt(numpy.sum(rs**2., axis=0))
    Ws = 0.5*v2s - mus/r
    a = -mus/2/Ws #semi-major axis

    L = numpy.array([rs[1,:]*vs[2,:] - rs[2,:]*vs[1,:],
                     rs[2,:]*vs[0,:] - rs[0,:]*vs[2,:],
                     rs[0,:]*vs[1,:] - rs[1,:]*vs[0,:]]) #angular momentum
    L2s = numpy.sum(L**2., axis=0)
    p = L2s/mus #semi-latus rectum
    e = numpy.sqrt(1 - p/a) #eccentricity

    #eccentric anomaly
    cosE = (1 - r/a)/e
    sinE = numpy.sum(rs*vs, axis=0)/(e*numpy.sqrt(mus*a))
    E = numpy.arctan2(sinE,cosE)

    #inclination
    sinI = numpy.sqrt(L[0,:]**2. + L[1,:]**2.)/numpy.sqrt(L2s)
    cosI = L[2,:]/numpy.sqrt(L2s)
    I = numpy.arctan2(sinI,cosI)

    #argument of pericenter
    sinwsinI = (((vs[0,:]*L[1,:] - vs[1,:]*L[0,:])/mus - 
                 rs[2,:]/r)/(e))
    coswsinI = ((numpy.sqrt(L2s)*vs[2,:])/mus -
                (L[0,:]*rs[1,:] - L[1,:]*rs[0,:])/(numpy.sqrt(L2s)*r))/(e)
    omega = numpy.arctan2(sinwsinI,coswsinI)

    #longitude of ascending node
    cosOsinI = -L[2,:]/(numpy.sqrt(L2s))
    sinOsinI = L[1,:]/(numpy.sqrt(L2s))
    Omega = numpy.arctan2(sinOsinI,cosOsinI)

    return a, e, E, I, omega, Omega


