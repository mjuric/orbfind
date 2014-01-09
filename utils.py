import math

# Input Units: vx, vy, vz (km s^-1) : ra, dec (deg)
# Output Units: km s^-1
def vel2pm(vx, vy, vz, ra, dec):
    c_ra = math.cos(ra * math.pi/180)
    s_ra = math.sin(ra * math.pi/180)
    c_dec = math.cos(dec * math.pi/180)
    s_dec = math.sin(dec * math.pi/180)

    v_ra  = -vx*s_ra  + vy*c_ra;
    tmp =  vx*c_ra  + vy*s_ra;
    vr  =  c_dec*tmp + vz*s_dec;
    v_dec  = -s_dec*tmp + vz*c_dec;

    return v_ra, v_dec, vr

# Input Units: km s^-1
# Output Units: km s^-1
def pm2vel(v_ra, v_dec, vr, ra, dec):
    c_ra = math.cos(ra * math.pi/180);
    s_ra = math.sin(ra *math.pi/180)
    c_dec = math.cos(dec * math.pi/180);
    s_dec = math.sin(dec * math.pi/180)

    tmp = s_dec*v_dec - c_dec*vr
    vx = -s_ra*v_ra -c_ra*tmp
    vy = c_ra*v_ra - s_ra*tmp
    vz = c_dec*v_dec + s_dec*vr

    return vx, vy, vz

def test():
    vx, vy, vz = pm2vel(0, 0, 10, 15, 25)
    v_ra, v_dec, vr = vel2pm(vx, vy, vz, 15, 25)
    print vx, vy, vz
    print v_ra, v_dec, vr

if __name__ == "__main__":
    test()
