from math import pi, sqrt, cos, factorial

import poliastro
from astropy import units as u
from astropy.units.quantity import Quantity
import numpy as np
import pandas as pd

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



def elliptical_period(a, meu):
    period = 2*pi*a**(1.5) / sqrt(meu)
    return period

#************************* Conservative Variables *****************************
def specific_energy_from_velo(velo, meu, radius):
    energy = ((velo**2)/2) - (meu/radius)
    return energy

def angular_momentum_from_p(p, meu):
    momentum = np.sqrt(p*meu) 
    return momentum

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
def time_of_flight_kepler(e, a, theta1, theta2, meu, pass_periapsis=0):
    E1 = get_eccentric_anomaly(e, theta1).value
    E2 = get_eccentric_anomaly(e, theta2).value
    n = np.sqrt(meu/a**3)

    dt = (2*np.pi*pass_periapsis + (E2 - e*np.sin(E2)) - (E1 - e*np.sin(E1))) / n
    return dt

def predict_location(e, a, theta1, dt, pass_periapsis, meu):
    n = np.sqrt(meu/a**3)
    E_o = get_eccentric_anomaly(e, theta1).value
    M_init = n.value*dt.value - 2*pass_periapsis*np.pi + (E_o - e*np.sin(E_o))
    guess_E = 0.5
    last_E = np.pi
    margin = 0.000001
    while abs(last_E - guess_E) > margin:
        M_current = guess_E - e*np.sin(guess_E)
        M_diff = M_init - M_current
        dMdE = 1 - e*np.cos(guess_E)
        dE = M_diff / dMdE
        last_E = guess_E
        guess_E = guess_E + dE

        print(f"En (rad): {last_E}")
        print(f"M-Mn: {M_diff}")
        print(f"dM/dE: {dMdE}")
        print(f"En+1: {guess_E}")
        print(abs(last_E - guess_E))
    theta2 = np.arccos((np.cos(guess_E) - e) / (1 - e*np.cos(guess_E)))
    theta2 = 2*np.pi - theta2
    return theta2*u.rad

def time_of_flight_universal_var(r_init, v_init, dt, meu, SandC=True):
    margin = 0.0001
    r0 = np.linalg.norm(r_init).to(DU_EARTH)
    v0 = np.linalg.norm(v_init).to(DUTU_EARTH)
    spec_energy = specific_energy_from_velo(v0, meu, r0).to(SPEC_E_EARTH)
    a = semi_major_axis_from_energy(spec_energy, meu).to(DU_EARTH)
    r_dot_v = np.dot(r_init, v_init).to(MOMENTUM_EARTH)

    x_guess = ((np.sqrt(meu)*dt)/a).to(DU_EARTH**(1/2))

    print(f"r0: {r0}")
    print(f"v0: {v0}")
    print(f"energy: {spec_energy}")
    print(f"a: {a}")
    print(f"r.v: {r_dot_v}")

    if SandC:
        if a > 0:
            SandC_func = get_SandC_elliptical
        elif a < 0: 
            SandC_func = get_SandC_hyperbolic
        if abs(get_z(x_guess, a)) < 0.001:
            SandC_func = get_SandC_parabolic

        counter = 0

        #while True:
        while counter < 5:
            print(f"x: {x_guess}")
            z = get_z(x_guess, a).value
            print(f"z: {z}")
            S, C = SandC_func(z)
            print(f"S,C: {S}, {C}")
            t = get_time_from_SandC(x_guess, z, S, C, r0, r_dot_v, meu) 
            print(f"t: {t}")           
            t_diff = dt - t
            print(f"tdiff: {t_diff}")
            if abs(t_diff).value < margin:
                x = x_guess
                break
            r = get_r_from_SandC(x_guess, z, S, C, r0, r_dot_v, meu)
            dtdx = r / np.sqrt(meu)
            x_guess = x_guess + t_diff/dtdx
            counter = counter + 1

    print(f"x_final: {x}")

    tof = 0
    r_final = 0
    v_final = 0
    return (tof, r_final, v_final)

# S and C variation Function
def get_SandC_elliptical(z):
    root_z = np.sqrt(z)
    S = (root_z - np.sin(root_z)) / np.sqrt(z**3)
    C = (1 - np.cos(root_z)) / z
    return (S, C)

def get_SandC_hyperbolic(z):
    neg_root_z = np.sqrt(-z)
    S = (np.sinh(neg_root_z) - neg_root_z) / np.sqrt((-z)**3)
    C = (1 - np.cosh(neg_root_z)) / z
    return (S, C)

def get_SandC_parabolic(z):
    min_stop = 0.000001
    exponent = 0
    sign = 1
    denom_S = 3
    denom_C = 2

    S = calc_z_series(z, exponent, sign, denom_S, min_stop)
    C = calc_z_series(z, exponent, sign, denom_C, min_stop)
    return (S, C)

def get_time_from_SandC(x, z, S, C, r_mag, r_dot_v, meu):
    term1 = x**3 * S
    term2 = (r_dot_v / np.sqrt(meu)) * x**2 * C
    term3 = r_mag * x * (1 - (z*S))
    time = (term1 + term2 + term3) / np.sqrt(meu)

    return time

def get_r_from_SandC(x, z, S, C, r0, r_dot_v, meu):
    term1 = x**2 * C
    term2 = (r_dot_v/np.sqrt(meu))*x*(1-(z*S))
    term3 = r0*(1 - (z*C))
    r = term1 + term2 + term3

    return r

def get_z(x, a):
    z = x**2 / a
    return z

# Universal Variable Normal Functions
def get_time_universal_var(r_mag, r_dot_v, meu, a, x):
    term1 = a*(x - np.sqrt(a)*np.sin(x/np.sqrt(a)))
    term2 = (r_dot_v/np.sqrt(meu))*a*(1 - np.cos(x/np.sqrt(a)))
    term3 = r_mag*np.sqrt(a)*np.sin(x/np.sqrt(a))
    t = (term1 + term2 + term3)/np.sqrt(meu)
    return t

# Misc
def get_eccentric_anomaly(e, theta):
    """ Returns Eccentric Anomaly
    theta (rad) - true anomaly """
    try:
        if theta.unit == "deg":
            theta = theta.to(u.rad)
    except:
        pass
    num = e + np.cos(theta)
    denom = 1 + e*np.cos(theta)
    E = np.arccos(num/denom)
    if theta.value > np.pi:
        E = 2*np.pi - E.value
    return E*u.rad

def calc_z_series(z, exponent, sign, denom, min):
    array = []
    hit_min = False
    while not hit_min:
        new_term = sign*(z**exponent)/factorial(denom)
        array.append(new_term)
        if new_term < min:
            hit_min = True
        sign = sign * -1
        denom = denom + 2
        exponent = exponent + 1
    sum = pd.Series(array).sum()
    return sum