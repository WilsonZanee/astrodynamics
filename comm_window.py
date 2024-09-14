import datetime as dt

import numpy as np
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.sampling import EpochsArray, TrueAnomalyBounds, EpochBounds
from poliastro.util import time_range
import satvis

a = 6778*u.km
ecc = 0.093315*u.one
inc = 1.85*u.deg
raan = 0*u.deg
argp = 0*u.deg
nu = 0*u.deg

start_date = Time("2026-01-01 00:00", scale="utc")
end_date = Time("2026-01-30 00:00", scale="utc")

min_per_month = 1440*30*u.min
t = np.arange(0, min_per_month.value)
print(t)

roadrunner = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)
rr_ephem = roadrunner.to_ephem(strategy=EpochsArray(epochs=time_range(start_date, end_date)))
roadrunner_r = []
for minute in t:




