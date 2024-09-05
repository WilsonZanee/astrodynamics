from math import pi, sqrt, cos

import poliastro
from astropy import units as u

MU_EARTH = u.def_unit("Mu Earth", 5.972e24*u.kg) 
DU_EARTH = u.def_unit("Du Earth", 6379.1*u.km) 
TU_EARTH = u.def_unit("Tu Earth", 806.8*u.s) 
DUTU_EARTH = u.def_unit("Du/Tu Earth", 7.9053661*u.km/u.s) 
MEU_EARTH = u.def_unit("meu Earth", 3.986e5*u.km**3/u.s**2)

AU_SUN = u.def_unit("Au Sun", 1.496e8*u.km) 
TU_SUN = u.def_unit("Tu Sun", 58.132821*u.d) 
AUTU_SUN = u.def_unit("Au/Tu Sun", 29.784852*u.km/u.s) 
MEU_SUN = u.def_unit("meu Sun", 1.3271544e11*u.km**3/u.s**2)


def elliptical_period(a, meu):
    period = 2*pi*a**(3/2) / sqrt(meu)
    return period

def specific_energy_from_velo(velo, meu, radius):
    energy = (velo**2/2) - (meu/radius)
    return energy

def angular_momentum_from_p(p, meu):
    momentum = sqrt(p*meu) 
    return momentum

def semi_major_axis_from_energy(energy, meu):
    a = -meu/(2*energy)
    return a

def eccentricity_from_momentum_energy(angular_momentum, energy, meu):
    e = sqrt(1 + (2*angular_momentum**2*energy/meu**2))
    return e

def orbit_radius_from_p_eccentricity_true_anomaly(e, p, theta):
    r = p/(1+e*cos(theta))
    return r
