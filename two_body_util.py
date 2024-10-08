from math import pi, sqrt, cos

import poliastro
from astropy import units as u
from astropy.units.quantity import Quantity
import numpy as np

MU_EARTH = u.def_unit("Mu Earth", 5.972e24*u.kg) 
DU_EARTH = u.def_unit("Du Earth", 6379.1*u.km) 
TU_EARTH = u.def_unit("Tu Earth", 806.8*u.s) 
DUTU_EARTH = u.def_unit("Du/Tu Earth", 7.9053661*u.km/u.s) 
MEU_EARTH = u.def_unit("meu Earth", 3.986e5*u.km**3/u.s**2)
SPEC_E_EARTH = u.def_unit("Du^2/Tu^2", DU_EARTH**2/TU_EARTH**2)
MOMENTUM_EARTH = u.def_unit("Du^2/Tu", DU_EARTH**2/TU_EARTH)

du_tu = (MEU_EARTH**(1/2)/DU_EARTH**(1/2), DUTU_EARTH, 
          lambda x: 1*x, lambda x: 1*x)

AU_SUN = u.def_unit("Au Sun", 1.496e8*u.km) 
TU_SUN = u.def_unit("Tu Sun", 58.132821*u.d) 
AUTU_SUN = u.def_unit("Au/Tu Sun", 29.784852*u.km/u.s) 
MEU_SUN = u.def_unit("meu Sun", 1.3271544e11*u.km**3/u.s**2)

spec_energy = u.km**2/u.s**2
angular_momentum = u.km**2/u.s

earth_radius_vector = np.matrix([[0],[0],[1]])*DU_EARTH
earth_rotational_velo = 7.2921159e-5*u.rad/u.s

earth_rotation_velo_vector = np.array([0, 0, earth_rotational_velo.value]
                                        )*earth_rotational_velo.unit

def elliptical_period(a, meu):
    period = 2*pi*a**(1.5) / sqrt(meu)
    return period

#************************* Conservative Variables *****************************
def specific_energy_from_velo(velo, meu, radius):
    energy = ((velo**2)/2) - (meu/radius)
    return energy

def specific_energy_from_rpra(rp, ra, meu):
    energy = -meu / (ra + rp)
    return energy.to(SPEC_E_EARTH)

def angular_momentum_from_p(p, meu):
    momentum = np.sqrt(p*meu) 
    return momentum

def angular_momentum_from_periapsis(vp, rp):
    h = vp*rp
    return h

#*********************** Velo Calcs *******************************************

def velo_from_energy(energy, meu, radius):
    velo = np.sqrt(2*(energy + meu/radius))
    return velo.to(DUTU_EARTH)

def velo_from_radius(meu, radius, semi_major_axis):
    velo = np.sqrt((2*meu/radius) - (meu/semi_major_axis))
    return velo.to(DUTU_EARTH)

def get_escape_velo(meu, radius):
    velo = np.sqrt(2*meu/radius)
    return velo.to(u.km/u.s)

def get_plane_change_dv(v1, v2, di):
    di = ensure_rad(di)
#    v1sqr = v1**2
#    v2sqr = v2**2
#    negterm = -2*v1*v2*cos(di.value)
#    dv = np.sqrt(v1sqr + v2sqr + negterm)
    dv = np.sqrt(v1**2 + v2**2 - 2*v1*v2*cos(di.value))
    return dv

def get_hyperbolic_excess_speed(velo_burnout, velo_esc):
    v_inf = np.sqrt(velo_burnout**2 - velo_esc**2)
    return v_inf

#************************ Orbital Elements ************************************
def semi_major_axis_from_energy(energy, meu):
    a = -meu/(2*energy)
    return a

def eccentricity_from_momentum_energy(angular_momentum, energy, meu):
    e = sqrt(1 + (2*angular_momentum**2*energy/meu**2))
    return e

#************************ Radius and Velocity Vectors *************************
def orbit_radius_from_p_eccentricity_true_anomaly(e, p, theta):
    r = p/(1+e*np.cos(theta))
    return r

#************************ Time of Flight *************************
def time_of_flight_kepler(e, meu, a, theta1, theta2, pass_periapsis=0):
    E1 = get_eccentric_anomaly(e, theta1)
    E2 = get_eccentric_anomaly(e, theta2)
    n = np.sqrt(meu/a**3)

    dt = (2*np.pi*pass_periapsis - (E2 - e*np.sin(E2)) - (E1 - e*np.sin(E1))) / n
    return dt


def get_eccentric_anomaly(e, theta):
    """ Returns Eccentric Anomaly
    theta (rad) - true anomaly """
    try:
        if isinstance(theta.units, u.deg):
            theta = theta.to(u.rad)
    except:
        pass
    num = e + np.cos(theta)
    denom = 1 + e*np.cos(theta)
    E = np.arccos(num/denom)
    if theta > np.pi:
        E = 2*np.pi - E
    return E

def ensure_rad(angle):
    if isinstance(angle, Quantity):
        if angle.unit == u.deg:
            rad_angle = angle.to(u.rad)
        elif angle.unit == u.rad:
            rad_angle = angle
        else:
            try:
                rad_angle = angle.to(u.rad)
            except:
                print(f"{angle} with unit {angle.unit} does not have "
                      f"a valid unit type")
    else:
        rad_angle = angle*u.rad
    return rad_angle



def get_sec(time_str):
    """Get seconds from time."""
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)