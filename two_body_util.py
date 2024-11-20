from math import pi, sqrt, cos, factorial, floor, isnan
from functools import partial

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
MEU_SUN = u.def_unit("mu Sun", 1.3271544e11*u.km**3/u.s**2)

LUNAR_RADIUS = 1737.4*u.km
LUNAR_MASS = 7.349e22*u.kg
MU_LUNAR = 4903*u.km**3/u.s**2
D_EARTH_LUNAR = 384400*u.km
LUNAR_SOI = 66300*u.km

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

def specific_energy_from_velo_infinity(velo):
    energy = ((velo**2)/2)
    return energy

def specific_energy_from_rpra(rp, ra, meu):
    energy = -meu / (ra + rp)
    return energy.to(u.km**2/u.s**2)

def angular_momentum_from_p(p, meu):
    momentum = np.sqrt(p*meu) 
    return momentum

def angular_momentum_from_periapsis(vp, rp):
    h = vp*rp
    return h

def angular_momentum_from_rv_angle(r, v, angle):
    h = r * v * np.cos(angle)
    return h

#*********************** Velo Calcs *******************************************

def velo_from_energy(energy, meu, radius):
    velo = np.sqrt(2*(energy + meu/radius))
    return velo.to(u.km/u.s)

def velo_from_radius(meu, radius, semi_major_axis):
    velo = np.sqrt((2*meu/radius) - (meu/semi_major_axis))
    return velo.to(u.km/u.s)

def get_escape_velo(meu, radius):
    velo = np.sqrt(2*meu/radius)
    return velo.to(u.km/u.s)

def get_plane_change_dv(v1, v2, di):
    di = ensure_rad(di)
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

def eccentricity_from_rarp(ra, rp):
    e = (ra - rp) / (ra + rp)
    return e

def p_from_h_mu(angular_momentum, mu):
    p = angular_momentum ** 2 / mu
    return p

def get_rp_from_p_e(p, e):
    rp = p / (1+e)
    return rp

#************************ Radius and Velocity Vectors *************************
def orbit_radius_from_p_eccentricity_true_anomaly(e, p, theta):
    r = p/(1+e*np.cos(theta))
    return r

#*********************** Interplanetary Trajectories **************************
def get_sphere_of_influence(large_mass, small_mass, distance):
    r = distance * (small_mass / large_mass) ** (2/5)
    return r

def get_flightpath_angle(angular_momentum, r, v, print_ang=False):
    cos_val = (angular_momentum.to(u.km**2/u.s) / (r.to(u.km) * v.to(u.km/u.s)))
    if cos_val > 1 or cos_val < -1:
        angle = 0*u.rad
    else:
        angle = np.arccos((angular_momentum.to(u.km**2/u.s) / (r * v)))
    if print_ang:
        print(f"h: {angular_momentum}")
        print(f"r,v: {r} {v}")
        print(f"cos_val: {cos_val}")
        print(f"Flightpath angle raw: {angle.to(u.deg)}")
    if angle > (np.pi / 2)*u.rad:
        angle = angle - (np.pi / 2)*u.rad
    elif angle < (-np.pi / 2)*u.rad:
        angle = angle + (np.pi / 2)*u.rad
    if print_ang:
        print(f"Flightpath angle final: {angle.to(u.deg)}")
    return angle

def get_v_inf(mu, r_planet, a_planet, transfer_angular_momentum, v_transfer,
              prints=False):
    v_planet = velo_from_radius(mu, r_planet, a_planet)
    flightpath_angle = get_flightpath_angle(transfer_angular_momentum, 
                                            r_planet, 
                                            v_transfer,
                                            print_ang=prints)
    v_inf = get_plane_change_dv(v_transfer, v_planet, flightpath_angle)
    if prints:
        print(f"v_planet: {v_planet}")
        print(f"v transfer: {v_transfer.to(u.km/u.s)}")
        print(f"v_inf: {v_inf.to(u.km/u.s)}")
    return v_inf

def get_best_turn_angle(rp, v_inf, mu, debug=False):
    if debug:
        print(rp)
        print(v_inf.to(u.km/u.s))
        print(mu)
    e = (1 + ((rp * v_inf**2) / mu))
    delta = 2*np.arcsin(1 / e)
    if debug:
        print(f"Max Delta: {delta.to(u.deg)}")
        print(f"e: {e}")
    return delta, e

def get_gravity_assist_velo(rp, v_transfer, v_planet, mu, debug=False):
    v_inf = v_transfer - v_planet
    v_inf_magnitude = np.linalg.norm(v_inf)
    d, e = get_best_turn_angle(rp, v_inf_magnitude, mu, debug=debug)
    R = np.array([[np.cos(d), -np.sin(d), 0],
                  [np.sin(d), np.cos(d), 0],
                  [0, 0, 1]])
    v_inf2 = np.dot(R, v_inf)
    v2_mag = np.linalg.norm(v_inf2)
    if debug:
        print(f"v2_mag: {v2_mag}")
        print(f"v2: {v_inf2}")
    v_final = v_inf2 + v_planet
    v_final_mag = np.linalg.norm(v_final)
    return v_final, v_final_mag

def get_offset_dist(rp, v_inf, mu_body):
    offset_dist = (rp / v_inf) * \
                    np.sqrt(v_inf ** 2 + (2 * mu_body / rp))
    return offset_dist

def angular_momentum_enter_SOI(r_target, v_inf, mu_body, delta=None):
    if delta is None:    
        delta = get_offset_dist(r_target, v_inf, mu_body)
    h = v_inf*delta
    return h

def get_SOI(dist, small_mass, big_mass):
    soi = dist * (small_mass / big_mass) ** (2/5)
    return soi

#**************************** Lunar Trajectories ******************************
def get_epsilon2(v_inf, r_inf, r_desired, mu):
    energy = specific_energy_from_velo(v_inf, mu, r_inf)
    vp = velo_from_energy(energy, mu, r_desired)
    h = angular_momentum_from_periapsis(vp, r_desired)
    phi = get_flightpath_angle(h, r_inf, v_inf)
    epsilon = (np.pi/2)*u.rad - phi
    return epsilon

def calc_epsilon2(v_moon, v_inf, vt_lunar, lambda1, phi1, gamma1):
    term1 = (v_moon / v_inf) * np.cos(lambda1)
    term2 = -(vt_lunar / v_inf) * np.cos(lambda1 - (phi1 - gamma1))
    epsilon2 = np.arcsin(term1 + term2)
    return epsilon2

def calc_rp_from_inf(epsilon, v_inf, r_inf, mu):
    energy = specific_energy_from_velo(v_inf, mu, r_inf)
    a = semi_major_axis_from_energy(energy, mu)
    if epsilon > 0:
        phi = np.pi/2*u.rad - epsilon
    else:
        phi = np.pi/2*u.rad + epsilon
    h = angular_momentum_from_rv_angle(r_inf, v_inf, phi)
    p = p_from_h_mu(h, mu)
    e = eccentricity_from_momentum_energy(h, energy, mu)
    rp = get_rp_from_p_e(p, e)
    return rp.to(u.km)

def get_radius_to_moon(lambda1):
    D = D_EARTH_LUNAR
    soi = LUNAR_SOI 
    r = np.sqrt(D**2 + soi**2 - 2*D*soi*np.cos(lambda1))
    return r

def get_gamma(soi, r1, lambda1):
    gamma = np.arcsin((soi / r1) * np.sin(lambda1))
    return gamma

def calc_rp_dv_from_patch_conic(lambda1, vt_lunar, r_earth_orbit, dv=False):
    r1 = get_radius_to_moon(lambda1)
    energy_transfer = specific_energy_from_velo(vt_lunar, 1*MEU_EARTH, r1)
    vt_earth = velo_from_energy(energy_transfer, 1*MEU_EARTH, r_earth_orbit)
    ht = angular_momentum_from_periapsis(vt_earth, r_earth_orbit)
    phi1 = get_flightpath_angle(ht, r1, vt_lunar)
    gamma1 = get_gamma(LUNAR_SOI, r1, lambda1)
    v_moon = velo_from_radius(1*MEU_EARTH, D_EARTH_LUNAR, D_EARTH_LUNAR)
    v_inf = get_plane_change_dv(vt_lunar, v_moon, phi1 - gamma1)
    epsilon2 = calc_epsilon2(v_moon, v_inf, vt_lunar, lambda1, phi1, gamma1)
    rp = calc_rp_from_inf(epsilon2, v_inf, LUNAR_SOI, MU_LUNAR)
    if dv:
        energy = specific_energy_from_velo(v_inf, MU_LUNAR, LUNAR_SOI)
        v_excess = velo_from_energy(energy, MU_LUNAR, rp)
        v_parked = velo_from_radius(MU_LUNAR, rp, rp)
        dv = v_excess - v_parked
        return(rp, epsilon2, dv)
    else:
        return (rp, epsilon2)


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
    #print(df) 
    theta2 = np.arccos((np.cos(guess_E*u.rad) - e) / (1 - e*np.cos(guess_E*u.rad)))
    if guess_E > np.pi and theta2.value < np.pi:
        theta2 = (2*np.pi - theta2.value)*u.rad
    elif guess_E < np.pi and theta2.value > np.pi:
        theta2 = (2*np.pi - theta2.value)*u.rad

#    if r_dot_v is not None:
#        if r_dot_v < 0:
#            theta2 = (2*np.pi - theta2.value)*u.rad
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
    #print(df) 
    f, g, f_dot, g_dot = get_fg(meu, x, z, S, C, r0, r, t)

    r_final = (f * r_init + g * v_init).to(r_init.unit)
    v_final = (f_dot * r_init + g_dot * v_init).to(v_init.unit)

    return (r_final, v_final)

# S and C variation Function
def get_SandC_elliptical(z):
    root_z = np.sqrt(z)
    try:
        if root_z.unit == "":
            root_z = root_z*u.rad
    except:
        root_z = root_z*u.rad
    S = (root_z.value - np.sin(root_z)) / np.sqrt(z**3)
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

    return time.to(u.day)

def get_r_from_SandC(x, z, S, C, r0, r_dot_v, meu):
    term1 = x**2 * C
    term2 = (r_dot_v/np.sqrt(meu))*x*(1-(z*S))
    term3 = r0*(1 - (z*C))
    r = term1 + term2 + term3

    return r

def get_fg(meu, x, z, S, C, r0, r, t):
    f = 1 - (x**2 / r0) * C
    g = (t - (x**3 / np.sqrt(meu)) * S).to(u.day)
    f_dot = (((np.sqrt(meu) * x) / (r0 * r)) * ((z * S) - 1)).to(1/u.day)
    g_dot = 1 - (x**2 / r) * C

    vars = {"f": f,
            "g": g,
            "f_dot": f_dot,
            "g_dot": g_dot}
    return (f, g, f_dot, g_dot)
    
def get_fg_gauss(r1, r2, y, x, z, S, A, meu):
    f = 1 - (y / r1).value
    g = A * np.sqrt(y / meu)
    f_dot = ((- np.sqrt(meu) * x) / (r1 * r2)) * (1 - (z * S))
    g_dot = 1 - (y / r2).value
    return (f, g, f_dot, g_dot)

def get_z(x, a):
    z = x**2 / a
    return z

# Gauss' Problem

def get_velo_gauss_problem(r1, r2, dt, meu, zguess=10, margin=1e-7,
                           max_iter=30, print_convergence_table=False):
    results = {}
    r1_0 = np.linalg.norm(r1)
    r2_0 = np.linalg.norm(r2)
    theta_short = np.arccos(np.dot(r1, r2) / (r1_0 * r2_0))
    theta_long = (2*np.pi)*u.rad - theta_short
    thetas = {"short": theta_short, "long": theta_long}

    z_guesses = np.linspace(1, 39, 6)
    z_guesses = np.insert(z_guesses, 0, zguess)

    for trip_type, d_theta in thetas.items():
        for guess in z_guesses:
            z, printout, converged = iteratively_calc_z_gauss(r1_0, 
                                                              r2_0,
                                                              d_theta,
                                                              guess,
                                                              dt,
                                                              meu,
                                                              margin,
                                                              max_iter)
            if converged:
                break
            else:
                if print_convergence_table:
                    print(f"{trip_type} did not converge with z guess {guess}")
                    print(printout)
        if print_convergence_table:
            print(f"{trip_type}:")
            print(printout)
        f, g, f_dot, g_dot = get_fg_gauss(r1_0, 
                                          r2_0, 
                                          printout["y"].iloc[-1],
                                          printout["x"].iloc[-1], 
                                          z, 
                                          printout["S"].iloc[-1], 
                                          printout["A"].iloc[-1],
                                          meu)

        v1 = ((r2 - (f * r1)) / g)
        v2 = ((g_dot * r2) - r1) / g

        results[trip_type] = {"v1": v1.value*(r1.unit/dt.unit),
                              "v2": v2.value*(r1.unit/dt.unit)}

    return results
            
def iteratively_calc_z_gauss(r1_0, r2_0, d_theta, zguess, dt, meu, margin, 
                             max_iter):

    x_list = []
    y_list = []
    z_list = []
    oldz_list = []
    S_list = []
    C_list = []
    A_list = []
    time_list = []
    dt_list = []
    dtdz_list = []
    
    if zguess > 0:
        SandC_func = get_SandC_elliptical
    elif zguess < 0: 
        SandC_func = get_SandC_hyperbolic
    if abs(zguess) < 1e-7:
        SandC_func = get_SandC_parabolic

    counter = 0
    z = zguess
    while counter < max_iter:
        A = get_A(r1_0, r2_0, d_theta)
        S, C = SandC_func(z)
        y, x = get_xy_gauss_problem(r1_0, r2_0, A, z, S, C)
        t = get_tof_gauss_problem(x, y, A, S, meu)
        t_diff = dt - t 
        dtdz = get_dtdz_gauss(z, S, C, A, y, x, meu)
        try:
            z_old = z.value
        except:
            z_old = z
        z = z + (t_diff / dtdz)

        if abs(t_diff).value < margin:
            break
        counter = counter + 1
        
        x_list.append(x.value)
        y_list.append(y.value)
        z_list.append(z.value)
        oldz_list.append(z_old)
        S_list.append(S)
        C_list.append(C)
        A_list.append(A.value)
        time_list.append(t.value)
        dt_list.append(t_diff.value)
        dtdz_list.append(dtdz.value)

    if abs(t_diff.value) > margin or isnan(z):
        converged = False
    else:
        converged = True

    printout = {
        "z_n": oldz_list,
        "x": x_list,
        "y": y_list,
        "S": S_list,
        "C": C_list,
        "A": A_list,
        "t": time_list,
        "dt": dt_list,
        "dtdz": dtdz_list,
        "z_n+1": z_list
    }
    printout = pd.DataFrame(printout)

    return (z, printout, converged) 


def get_tof_gauss_problem(x, y, A, S, meu):
    num = (x**3 * S) + (A * np.sqrt(y))
    denom = np.sqrt(meu)
    tof = num / denom
    return tof
    
def get_xy_gauss_problem(r1_0, r2_0, A, z, S, C):
    num = 1 - (z * S)
    denom = np.sqrt(C)
    y = r1_0 + r2_0 - A * num / denom
    x = np.sqrt(y/C)

    return (y, x)

def get_A(r1_mag, r2_mag, d_theta):
    num = np.sqrt(r1_mag * r2_mag) * np.sin(d_theta)
    denom = np.sqrt(1 - np.cos(d_theta))
    A = num / denom
    return A

def get_dtdz_gauss(z, S, C, A, y, x, meu):
    dSdz = get_dSdz(z, S, C)
    dCdz = get_dCdz(z, S, C)

    term1 = x**3 * (dSdz - ((3 * S * dCdz) / (2 * C)))
    term2 = (A / 8) * (((3 * S * np.sqrt(y)) / C) + (A / x))
    dtdz = (term1 + term2) / np.sqrt(meu)
    return dtdz

def get_dSdz(z, S, C):
    dSdz = (C - (3 * S)) / (2 * z)
    return dSdz

def get_dCdz(z, S, C):
    dCdz = (1 - (z * S) - 2 * C) / (2 * z)
    return dCdz

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

def bisection_method(lower_bound, upper_bound, function, desire_result, 
                     accuracy=1e-6, debug=False):
    y_mid = 1
    while abs(y_mid) > accuracy:
        y_lower = function(lower_bound) - desire_result
        y_upper = function(upper_bound) - desire_result
        if y_lower > 0 and y_upper < 0:
            up = upper_bound
            upper_bound = lower_bound
            lower_bound = up
        midpt = (lower_bound + upper_bound) / 2
        y_mid = function(midpt) - desire_result
        if debug:
            print(f"low: {lower_bound}, yl: {y_lower}")
            print(f"upper: {upper_bound}, yu: {y_upper}")
            print(f"midpt: {midpt}, y_mid: {y_mid}")
        if y_mid > 0:
            upper_bound = midpt
        else:
            lower_bound = midpt
    return midpt


