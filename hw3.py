import os

import numpy as np
from astropy import units as u
from astropy.units import imperial as imp

from orbits import Orbit, OrbitalElements, RadarObservation

import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: Radar observation to Orbital Elements ---------------------------
print("\n" + "-"*width + "\n")
print("Question 1: Radar observation to Orbital Elements\n")

#range = 3*util.DU_EARTH
#range_rate = 0.4*util.DUTU_EARTH
#az = (30*u.deg).to(u.rad)
#az_rate = 0*u.rad/util.TU_EARTH
#el = (10*u.deg).to(u.rad)
#el_rate = 0.1*u.rad/util.TU_EARTH

range = 0.4*util.DU_EARTH
range_rate = 0*util.DUTU_EARTH
az = (90*u.deg).to(u.rad)
az_rate = 10*u.rad/util.TU_EARTH
el = (30*u.deg).to(u.rad)
el_rate = 5*u.rad/util.TU_EARTH

ref_zulu_sid_time = 0*u.deg
lat = 60*u.deg
long = -150*u.deg
dt = util.get_sec("06:00:00")*u.s


obs = RadarObservation(range, 
                       range_rate, 
                       az, 
                       az_rate, 
                       el, 
                       el_rate, 
                       lat=lat, 
                       long=long,
                       ref_zulu_theta=ref_zulu_sid_time, 
                       dt=dt)

r, v = obs.get_vectors_in_geocentric_IJK()
print(f"r{r}")
print(f"v{v}")

meu = 1*util.MEU_EARTH

orbit = Orbit(r_vector=r, v_vector=v, meu=meu)
print(orbit.orbital_elements)

# Question 2: Orbit Circularization------------------------------------
print("\n" + "-"*width + "\n")
print("Question 2: Orbit Circularization\n")

meu = 1*util.MEU_EARTH
rp = (1*util.DU_EARTH + 100*u.km)
ra = (1*util.DU_EARTH + 250*u.km)
energy1 = util.specific_energy_from_rpra(ra, rp, meu)
velo1 = util.velo_from_energy(energy1, meu, ra)
velo2 = util.velo_from_radius(meu, ra, ra)
dv = velo2 - velo1

print(f"Delta V Required: {dv} = {dv.to(u.km/u.s)}")

# Question 3: Escape to Mars -------------------------------------------------
print("\n" + "-"*width + "\n")
print("Question 3: Escape to Mars\n")

meu = 1*util.MEU_EARTH
circular_radius = (1*util.DU_EARTH + 200*u.km)
dv = 12*u.km/u.s

parked_orbit_velo = util.velo_from_radius(meu, circular_radius, circular_radius)
print(f"a. Parked ORbit Velo: {parked_orbit_velo} ="
      f" {parked_orbit_velo.to(u.km/u.s)}")
v_burnout = parked_orbit_velo + dv
print(v_burnout)
v_esc = util.get_escape_velo(meu, circular_radius)
v_inf = util.get_hyperbolic_excess_speed(v_burnout, v_esc)
print(v_inf)
