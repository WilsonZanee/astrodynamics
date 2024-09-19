import numpy as np
from numpy import cos, sin

import two_body_util as util

IJK = "ijk" # Geocentric
SEZ = "sez" # Topocentric Horizon
PQW = "pqw" # Parifocal

def perifocal_to_geocentric_matrix(raan, arg_o_periapsis, inclination):
    omega = arg_o_periapsis
    i = inclination
    cord11 = cos(raan)*cos(omega) - sin(raan)*sin(omega)*cos(i)
    cord12 = -cos(raan)*sin(omega) - sin(raan)*cos(omega)*cos(i)
    cord13 = sin(raan)*sin(i)
    cord21 = sin(raan)*cos(omega) + cos(raan)*sin(omega)*cos(i)
    cord22 = -sin(raan)*sin(omega) + cos(raan)*cos(omega)*cos(i)
    cord23 = -cos(raan)*sin(i)
    cord31 = sin(omega)*sin(i)
    cord32 = cos(omega)*sin(i)
    cord33 = cos(i)
    r = np.matrix([[cord11, cord12, cord13],
                   [cord21, cord22, cord23],
                   [cord31, cord32, cord33]])
    return r

def topocentric_to_geocentrix_matrix(lat, local_siderial_t):
    lat = util.ensure_rad(lat)
    theta = util.ensure_rad(local_siderial_t)
    D = np.matrix([[sin(lat)*cos(theta), -sin(theta), cos(lat)*cos(theta)],
                   [sin(lat)*sin(theta), cos(theta), cos(lat)*sin(theta)],
                   [-cos(lat), 0, sin(lat)]])
    return D

def get_local_siderail_time(ref_zulu_siderial_time, dt, lat):
    theta_go = util.ensure_rad(ref_zulu_siderial_time)
    lat = util.ensure_rad(lat)
    local_sid_time = (theta_go 
                      + util.earth_rotational_velo*dt 
                      + lat)
    return local_sid_time