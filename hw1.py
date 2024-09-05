import datetime as t

from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.constants.general import GM_earth, GM_sun, R_mean_earth
from astropy import units as u
from astropy import constants as C

import two_body_util as util

# Problem 2: Halley's Commet

period = util.elliptical_period(17.9564, 1)
print(f"Orbital Period in canonical units: {period}")
orbit_start_time = t.datetime(1986, 2, 9)
period_days = period*util.TU_SUN.to(u.d)
time_delta = t.timedelta(days=period_days)
orbit_end_time = orbit_start_time + time_delta

print(f"Orbital Period in days: {period_days}")
print(f"Halley's comet will return to perihelion: {orbit_end_time}")

# Question 3: Earth Satellite

meu = 1*util.MEU_EARTH
print(meu)
p = 13778*(u.km)
p = p.to(util.DU_EARTH)
print(p)
velo_p = 6.5*(u.km/u.s)
velo_p = velo_p.to(util.DUTU_EARTH)
print(velo_p)
momentum = util.angular_momentum_from_p(p.value, meu.value)
print(f"Momentum: {momentum}")
spec_energy = util.specific_energy_from_velo(velo_p.value, meu.value, p.value)
print(f"Specific Energy in Canonical Units: {spec_energy}")
a = util.semi_major_axis_from_energy(spec_energy, meu.value)
print(f"a: {a}")
e = util.eccentricity_from_momentum_energy(momentum, spec_energy, meu.value)
print(f"Eccentricity: {e}")



