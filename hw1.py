import datetime as t

from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.constants.general import GM_earth, GM_sun, R_mean_earth
from astropy import units as u
from astropy import constants as C

#MEU_SUN = 


import two_body_util as util

# Problem 2: Halley's Commet

period = util.elliptical_period(17.9564, 1)
print(f"Orbital Period in canonical units: {period}")
orbit_start_time = t.datetime(1986, 2, 9)
period_days = period * util.TU_SUN
time_delta = t.timedelta(days=period_days)
orbit_end_time = orbit_start_time + time_delta

print(f"Orbital Period in days: {period_days}")
print(f"Halley's comet will return to perihelion: {orbit_end_time}")

# Question 3: Earth Satellite

km_per_s = u.km/u.s
Du = u.def_unit("Du", R_mean_earth.to(u.km))
Tu = u.def_unit("Tu", 806.8*u.s)
meu_earth = 1

p = 13778*u.km.to(Du)
velo_p = 6.5*km_per_s #km/s
velo_p_canon = velo_p.to(Du/Tu)
momentum = util.angular_momentum_from_p(p, meu_earth)
print(momentum)
spec_energy = util.specific_energy_from_velo(velo_p_canon, 1, p)
print(spec_energy)


