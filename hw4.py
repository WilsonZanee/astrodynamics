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

orbit = Orbit(r_vector=r2, v_vector=v2, meu=meu_earth)
oe = orbit.orbital_elements
print(oe)

theta_final = util.predict_location(oe.e, oe.a, oe.theta, dt2, 1, meu_earth)

print(f"{theta_final}")

# Question 3: Kepler's Problem Solved w/ by S and C Method --------------------
print("\n" + "-"*width + "\n")
print("Question 3: Kepler's Problem Solved w/ by S and C Method\n")