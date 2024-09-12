import numpy as np
from numpy import cos, sin

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
    cord23 = -cos(raan)*cos(i)
    cord31 = sin(omega)*sin(i)
    cord32 = cos(omega)*sin(i)
    cord33 = cos(i)
    r = np.matrix([[cord11, cord12, cord13],
                   [cord21, cord22, cord23],
                   [cord31, cord32, cord33]])
    return r