from math import pi, sqrt, cos, factorial, floor

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

earth_radius_vector = np.matrix([[0],[0],[1]])*DU_EARTH
earth_rotational_velo = 7.2921159e-5*u.rad/u.s

earth_rotation_velo_vector = np.array([0, 0, earth_rotational_velo.value]
                                        )*earth_rotational_velo.unit

def elliptical_period(a, meu):
    period = 2*pi*a**(1.5) / np.sqrt(meu)
    return period

def get_pass_periapsis(e, a, theta, tof, meu):
    theta_final = 2*np.pi*u.rad
    dt = time_of_flight_kepler(e, a, theta,theta_final, meu).to(u.s)
    period = elliptical_period(a, meu).to(u.s)

    k = floor((dt - tof.to(u.s)) / period)
    return k

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
def time_of_flight_kepler(e, a, theta1, theta2, meu, pass_periapsis=0):
    E1 = get_eccentric_anomaly(e, theta1).value
    E2 = get_eccentric_anomaly(e, theta2).value
    n = np.sqrt(meu/(a.to(DU_EARTH))**3).to(1/u.s)

    dt = (2*np.pi*pass_periapsis + (E2 - e*np.sin(E2)) - (E1 - e*np.sin(E1))) / n
    return dt

def predict_location(e, a, theta1, dt, pass_periapsis, 
                     meu, guess_E=2, r_dot_v=-1):
    n = np.sqrt(meu/a**3)
    E_o = get_eccentric_anomaly(e, theta1).value
    M_init = n*dt - 2*pass_periapsis*np.pi + (E_o - e*np.sin(E_o))
    last_E = np.pi
    margin = 0.000001
    e_list = []
    M_diff_list = []
    dMdE_list = []
    eNew_list = []
    while abs(last_E - guess_E) > margin:
        M_current = guess_E - e*np.sin(guess_E*u.rad)
        M_diff = M_init - M_current
        dMdE = 1 - e*np.cos(guess_E*u.rad)
        dE = M_diff / dMdE
        last_E = guess_E
        guess_E = guess_E + dE

        e_list.append(last_E)
        M_diff_list.append(M_diff)
        dMdE_list.append(dMdE)
        eNew_list.append(guess_E)
    
    printout = {
    "En": e_list,
    "M-Mn": M_diff_list,
    "dM/dE": dMdE_list,
    "En+1": eNew_list
    }
    df = pd.DataFrame(printout)
    print(df) 
    theta2 = np.arccos((np.cos(guess_E*u.rad) - e) / (1 - e*np.cos(guess_E*u.rad)))
    if r_dot_v is not None:
        if r_dot_v < 0:
            theta2 = (2*np.pi - theta2.value)*u.rad
    return theta2

def time_of_flight_universal_var(r_init, v_init, dt, meu, SandC=True, 
                                 max_iter=30, hyperbolic_guess=2):
    margin = 1e-7
    r0 = np.linalg.norm(r_init).to(DU_EARTH)
    v0 = np.linalg.norm(v_init).to(DUTU_EARTH)
    spec_energy = specific_energy_from_velo(v0, meu, r0).to(SPEC_E_EARTH)
    a = semi_major_axis_from_energy(spec_energy, meu).to(DU_EARTH)
    r_dot_v = np.dot(r_init, v_init).to(MOMENTUM_EARTH)
    if a > 0:
        x_guess = ((np.sqrt(meu)*dt)/a).to(DU_EARTH**(1/2))
    if a < 0:
        x_guess = hyperbolic_guess*DU_EARTH**(1/2)

    x_list = []
    z_list = []
    S_list = []
    C_list = []
    time_list = []
    dt_list = []
    dtdx_list = []

    if SandC:
        if a > 0:
            SandC_func = get_SandC_elliptical
        elif a < 0: 
            SandC_func = get_SandC_hyperbolic
        if abs(get_z(x_guess, a)) < 1e-7:
            SandC_func = get_SandC_parabolic
        counter = 0

        #while True:
        while counter < max_iter:
            z = get_z(x_guess, a).value
            S, C = SandC_func(z)
            t = get_time_from_SandC(x_guess, z, S, C, r0, r_dot_v, meu) 
            t_diff = dt - t
            r = get_r_from_SandC(x_guess, z, S, C, r0, r_dot_v, meu)
            dtdx = r / np.sqrt(meu)
            old_x = x_guess
            
            x_list.append(old_x.value)
            z_list.append(z)
            S_list.append(S)
            C_list.append(C)
            time_list.append(t.value)
            dt_list.append(t_diff.value)
            dtdx_list.append(dtdx.value)    

            if abs(t_diff).value < margin:
                x = x_guess
                break
            x_guess = x_guess + t_diff/dtdx
            counter = counter + 1
    printout = {
    "x": x_list,
#    "z": z_list,
#    "S": S_list,
#    "C": C_list,
#    "time": time_list,
    "dt": dt_list,
    "dt/dx": dtdx_list
    }
    df = pd.DataFrame(printout)
    print(df) 
    f, g, f_dot, g_dot = get_fg(meu, x, z, S, C, r0, r, t)

    r_final = f * r_init + g * v_init
    v_final = f_dot * r_init + g_dot * v_init

    return (r_final, v_final)

# S and C variation Function
def get_SandC_elliptical(z):
    root_z = np.sqrt(z)
    S = (root_z - np.sin(root_z)) / np.sqrt(z**3)
    C = (1 - np.cos(root_z)) / z
    return (S, C)

def get_SandC_hyperbolic(z):
    z = np.longdouble(z)
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

def get_time_from_SandC(x, z, S, C, r0, r_dot_v, meu):
    term1 = (x**3) * S
    term2 = (r_dot_v / np.sqrt(meu)) * (x**2) * C
    term3 = r0 * x * (1 - (z*S))
    time = (term1 + term2 + term3) / np.sqrt(meu)

    return time.to(TU_EARTH)

def get_r_from_SandC(x, z, S, C, r0, r_dot_v, meu):
    term1 = x**2 * C
    term2 = (r_dot_v/np.sqrt(meu))*x*(1-(z*S))
    term3 = r0*(1 - (z*C))
    r = term1 + term2 + term3

    return r

def get_fg(meu, x, z, S, C, r0, r, t):
    f = 1 - (x**2 / r0) * C
    g = (t - (x**3 / np.sqrt(meu)) * S).to(TU_EARTH)
    f_dot = (((np.sqrt(meu) * x) / (r0 * r)) * ((z * S) - 1)).to(1/TU_EARTH)
    g_dot = 1 - (x**2 / r) * C

    vars = {"f": f,
            "g": g,
            "f_dot": f_dot,
            "g_dot": g_dot}

    return (f, g, f_dot, g_dot)

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
