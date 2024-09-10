import datetime as t
import os

import numpy as np
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.constants.general import GM_earth, GM_sun, R_mean_earth
from astropy import units as u
from astropy import constants as C

import two_body_util as util
from orbits import Orbit

width = os.get_terminal_size().columns

# Problem 2: Halley's Commet --------------------------------------------------
print("\n" + "-"*width)
print("Problem 2: Halley's Commet\n")


period = util.elliptical_period(17.9564, 1)
print(f"Orbital Period in canonical units: {period*util.TU_SUN}")
orbit_start_time = t.datetime(1986, 2, 9)
period_days = (period*util.TU_SUN).to(u.d)
time_delta = t.timedelta(days=period_days.value)
orbit_end_time = orbit_start_time + time_delta

print(f"Orbital Period in days: {period_days}")
print(f"Halley's comet will return to perihelion: {orbit_end_time}")


# Question 3: Earth Satellite -------------------------------------------------
print("\n" + "-"*width)
print("Question 3: Earth Satellite\n")

meu = (1*util.MEU_EARTH).to(u.km**3/u.s**2)
p = (13778*u.km)#.to(util.DU_EARTH)
velo_p = (6.5*(u.km/u.s))#.to(util.DUTU_EARTH)
momentum = util.angular_momentum_from_p(p, meu)
print(f"Momentum: {momentum} or {momentum.to(util.MOMENTUM_EARTH)}")
spec_energy = util.specific_energy_from_velo(velo_p, meu, p)
print(f"Specific Energy in Canonical Units: {spec_energy} or {spec_energy.to(util.SPEC_E_EARTH)}")
a = util.semi_major_axis_from_energy(spec_energy, meu)
print(f"a: {a}")
e = util.eccentricity_from_momentum_energy(momentum, spec_energy, meu)
print(f"Eccentricity: {e}")
theta = (100*u.deg).to(u.rad)
r = util.orbit_radius_from_p_eccentricity_true_anomaly(e, p, theta.value)
print(f"Radius of orbit at true anomaly {np.rad2deg(theta)*u.deg} is {r} or {r.to(util.DU_EARTH)}")


# Question 4: Satellite OE ----------------------------------------------------
print("\n" + "-"*width)
print("Question 4: Satellite OE\n")

r = np.array([6578, 6778, 6178])*u.km
v = np.array([-4.5, 2.5, 5.5])*u.km/u.s
meu = (1*util.MEU_EARTH).to(u.km**3/u.s**2)

orbit = Orbit(r_vector=r, v_vector=v, meu=meu)
print(f"Angular Momentum Vector: {orbit.angular_momentum}"
      f" = {orbit.angular_momentum.to(util.MOMENTUM_EARTH)}")
print(f"Specific Energy: {orbit.spec_energy} or {orbit.spec_energy.to(util.SPEC_E_EARTH)}")
print(f"Eccentricity: {orbit.eccentricity}")
print(f"rp: {orbit.rp} or {orbit.rp.to(util.DU_EARTH)}")
print(f"ra: {orbit.ra} or {orbit.ra.to(util.DU_EARTH)}")



