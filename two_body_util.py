from math import pi, sqrt

import poliastro

MU_EARTH = 5.972e24 # kg
DU_EARTH = 6379.1 # km
TU_EARTH = 806.8 # seconds
MEU_EARTH = 3.986e5 # kg^3/s^2

TU_SUN = 58.13 # Days

def elliptical_period(a, meu):
    period = 2*pi*a**(3/2) / sqrt(meu)
    return period

def specific_energy_from_velo(velo, meu, radius):
    energy = velo^2/2 - meu/radius
    return energy

def angular_momentum_from_p(p, meu):
    momentum = sqrt(p*meu) 
    return momentum


#def eccentricity():


