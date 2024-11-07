import os

import numpy as np
from astropy import units as u
from astropy.units import imperial as imp

from orbits import Orbit, OrbitalElements, RadarObservation

import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: Earth-Jupiter interplanetary transfer ---------------------------
print("\n" + "-"*width + "\n")
print("Question 1: Earth-Jupiter interplanetary transfer\n")

e_transfer = 0.7
r_earth = 1*util.AU_SUN
r_jupiter = 5.2*util.AU_SUN
mu_sun = 1*util.MEU_SUN

ra = Orbit.get_ra(r_earth, e_transfer)

# Question 1a: 
print(f"The transfer orbit with an eccentricity of 0.7 should work because\n "
      "assuming that earth is at the periapsis, the apoapsis of the oribit\n "
      f"is {ra} which is more than jupiters orbit radius of {r_jupiter}")

# Question 2a:

# Transfer Phase
print((r_earth + ra).to(u.km))
transfer_energy = util.specific_energy_from_rpra(r_earth, ra, mu_sun)\
                    .to(u.km**2 / u.s**2)
print(transfer_energy)
vt1 = util.velo_from_energy(transfer_energy, mu_sun, r_earth)
vt2 = util.velo_from_energy(transfer_energy, mu_sun, r_jupiter)

# Earth Escape Phase
v_earth = util.velo_from_radius(mu_sun, r_earth, r_earth)
v_inf1 = vt1 - v_earth

#excape_energy = util.specific_energy_from_velo()

print(v_inf1)

mu_earth = (1*util.MEU_EARTH).to(u.km**3 / u.s**2)
earth_orbit_r = 200*u.km + 1*util.DU_EARTH
v_earth_orbit = (util.velo_from_radius(mu_earth, earth_orbit_r, earth_orbit_r)\
                 .to(u.km/u.s))
mu_jupiter = 1.26713e8 * (u.km**3 / u.s**2)