import os

import numpy as np
from astropy import units as u
from astropy.units import imperial as imp

from orbits import Orbit, OrbitalElements, RadarObservation

import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: S/C in Saturn orbit TOF Calculation -----------------------------
print("\n" + "-"*width + "\n")
print("Question 1: S/C in Saturn orbit TOF Calculation\n")

e = 0.6
a = 9e5 * u.km
meu_saturn = (3.7931187e16 * u.m**3 / u.s**2).to(u.km**3 / u.s**2)

theta0 = (300*u.deg).to(u.rad)
theta1 = (190*u.deg).to(u.rad)
k = 1

dt1 = util.time_of_flight_kepler(e, a, theta0, theta1, meu_saturn, 
                                pass_periapsis=k)

print(f"Time of Flight: {dt1.to(u.hr)}")


# Question 2: S/C in Earth orbit True Anomaly from TOF Calculation ------------
print("\n" + "-"*width + "\n")
print("Question 2: S/C in Earth orbit True Anomaly from TOF Calculation\n")

r2 = [2, -1.287, -0.3]*util.DU_EARTH
v2 = [0.3, -0.63, 0.229]*util.DUTU_EARTH
dt2 = (2.22*u.hr).to(util.TU_EARTH)
meu_earth = 1*util.MEU_EARTH
r_dot_v = np.dot(r2, v2)

orbit = Orbit(r_vector=r2, v_vector=v2, meu=meu_earth)
oe = orbit.orbital_elements
k2 = util.get_pass_periapsis(oe.e, oe.a, oe.theta, dt2, meu_earth)

theta_final = util.predict_location(oe.e, oe.a, oe.theta, dt2, k2, meu_earth,
                                     r_dot_v=r_dot_v)
print(f"The final true anomaly is {theta_final} or {theta_final.to(u.deg)}")

print(util.time_of_flight_kepler(oe.e, oe.a, oe.theta, theta_final, meu_earth, pass_periapsis=k2).to(u.hr))


# Question 3: Kepler's Problem Solved w/ by S and C Method --------------------
print("\n" + "-"*width + "\n")
print("Question 3: Kepler's Problem Solved w/ by S and C Method\n")

# 3a

r3a = [0, 1, 0]*util.DU_EARTH
v3a = [0, 0, 1]*util.DUTU_EARTH
dt3a = 3.1315926*util.TU_EARTH

print("3a:")
r3a_f, v3a_f= util.time_of_flight_universal_var(r3a, v3a, dt3a, meu_earth)
print(f"\nr: {r3a_f}")
print(f"v: {v3a_f}")

# 3b

r3b = [0.5, 0.7, 0.8]*util.DU_EARTH
v3b = [0, 0.1, 0.9]*util.DUTU_EARTH
dt3b = -20*util.TU_EARTH

print("3b:")
r3b_f, v3b_f= util.time_of_flight_universal_var(r3b, v3b, dt3b, meu_earth)
print(f"\nr: {r3b_f}")
print(f"v: {v3b_f}")

# 3c

r3c = [0, 0, -0.5]*util.DU_EARTH
v3c = [0, 2, 0]*util.DUTU_EARTH
dt3c = (10e6)*util.TU_EARTH

print("3c:")
try:
    r3c_f, v3c_f= util.time_of_flight_universal_var(r3c, v3c, dt3c, meu_earth)
    print(f"\nr: {r3c_f}")
    print(f"v: {v3c_f}")
except Exception as e:
    print(e)
# 3d

r3d = [-0.5, 0, 0]*util.DU_EARTH
v3d = [0, 1.999, 0]*util.DUTU_EARTH
dt3d = 10e3*util.TU_EARTH

print("3d:")
r3d_f, v3d_f= util.time_of_flight_universal_var(r3d, v3d, dt3d, meu_earth)
print(f"\nr: {r3d_f}")
print(f"v: {v3d_f}")