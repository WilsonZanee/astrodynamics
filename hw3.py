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

range = 3*util.DU_EARTH
range_rate = 0.4*util.DUTU_EARTH
az = (30*u.deg).to(u.rad)
az_rate = 0*u.rad/util.TU_EARTH
el = (10*u.deg).to(u.rad)
el_rate = 0.1*u.rad/util.TU_EARTH

ref_zulu_sid_time = 0*u.deg
lat = 32.5*u.deg
long = -106.6*u.deg
dt = util.get_sec("18:00:00")*u.s

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

# Question 3: Escape to Mars -------------------------------------------_------
print("\n" + "-"*width + "\n")
print("Question 3: Escape to Mars\n")

# Part a
meu = 1*util.MEU_EARTH
circular_radius = (1*util.DU_EARTH + 200*u.km)

parked_orbit_velo = util.velo_from_radius(meu, circular_radius, circular_radius)
print(f"a. Parked Orbit Velo: {parked_orbit_velo} ="
      f" {parked_orbit_velo.to(u.km/u.s)}")

# Part b
v_burnout = (12*u.km/u.s).to(util.DU_EARTH/util.TU_EARTH)
energy = util.specific_energy_from_velo(v_burnout, meu, circular_radius)
a = util.semi_major_axis_from_energy(energy, meu)
angular_momentum = util.angular_momentum_from_periapsis(v_burnout, 
                                                        circular_radius)
eccentricity = util.eccentricity_from_momentum_energy(angular_momentum, 
                                                      energy, 
                                                      meu)

print(f"b. ")
print(f"    Eccentricity: {eccentricity}")
print(f"    Specific Energy: {energy.to(util.SPEC_E_EARTH)}")
print(f"    Angular Momentum: {angular_momentum.to(util.MOMENTUM_EARTH)}")

# Part c
v_esc = util.get_escape_velo(meu, circular_radius)
v_inf = util.get_hyperbolic_excess_speed(v_burnout, v_esc)
print(f"c. Hyperbolic Excess Speed: {v_inf} = {v_inf.to(u.km/u.s)}")

# Part d
velo_at_large_r = util.velo_from_energy(energy, meu, 1e6*u.km)
print(f"d. Velo at 1*10^6 km on Hyperbolic trajectory: {velo_at_large_r}")

# Question 4: Orbital Maneuver Comparison -------------------------------------
print("\n" + "-"*width + "\n")
print("Question 4: Orbital Maneuver Comparison\n")

meu = 1*util.MEU_EARTH
i = 28.5*u.deg
r_LEO = 1*util.DU_EARTH + 200*u.km
v_LEO = util.velo_from_radius(meu, r_LEO, r_LEO)
energy_LEO = util.specific_energy_from_velo(v_LEO, meu, r_LEO)
r_GEO = 1*util.DU_EARTH + 35786*u.km
v_GEO = util.velo_from_radius(meu, r_GEO, r_GEO)
energy_GEO = util.specific_energy_from_velo(v_GEO, meu, r_GEO)

rp_h = r_LEO
ra_h = r_GEO
energy_h = util.specific_energy_from_rpra(rp_h, ra_h, meu)
v_t1 = util.velo_from_energy(energy_h, meu, rp_h)
v_t2 = util.velo_from_energy(energy_h, meu, ra_h)

# Part a
# Circular LEO to hohmann elliptical
dv_a = v_t1 - v_LEO
# elliptical hohmann to circular GEO
dv_a = dv_a + (v_GEO - v_t2)
# Plane change to equatorial
dv_a = dv_a + util.get_plane_change_dv(v_GEO, v_GEO, i)
print(f"a. dv = {dv_a.to(u.km/u.s)}")

# Part b
# Circular LEO to hohmann elliptical
dv_a = v_t1 - v_LEO
# Plane change to equatorial
dv_a = dv_a + util.get_plane_change_dv(v_t2, v_t2, i)
# elliptical hohmann to circular GEO
dv_a = dv_a + (v_GEO - v_t2)
print(f"b. dv = {dv_a.to(u.km/u.s)}")

# Part c
# Circular LEO to hohmann elliptical
dv_a = v_t1 - v_LEO
# Plane change and burn to circular GEO
dv_a = dv_a + util.get_plane_change_dv(v_t2, v_GEO, i)
print(f"c. dv = {dv_a.to(u.km/u.s)}")
