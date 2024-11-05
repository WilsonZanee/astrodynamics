import os

import numpy as np
import astropy.units as u

from orbits import Orbit, OrbitalElements
import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: Gauss' Problem -----------------------------
print("\n" + "-"*width + "\n")
print("Question 1: Gauss' Problem\n")
print("1a")
r1 = [0.3, 0.7, 0.4]*util.DU_EARTH
r2 = [0.6, -1.4, 0.8]*util.DU_EARTH
dt = 5*util.TU_EARTH
meu = 1*util.MEU_EARTH

results = util.get_velo_gauss_problem(r1, r2, dt, meu, zguess=5, 
                                      print_convergence_table=True)
print(results["short"])

print("1b")
r1 = [0.5, 0.6, 0.7]*util.DU_EARTH
r2 = [0, 1, 0]*util.DU_EARTH
dt = 1.2*util.TU_EARTH
meu = 1*util.MEU_EARTH

results = util.get_velo_gauss_problem(r1, r2, dt, meu, zguess=5, 
                                      print_convergence_table=True )
print(results["long"])

# Question 2: Interplanetary Transfer -----------------------------
print("\n" + "-"*width + "\n")
print("Question 2: Interplanetary Transfer\n")

r_earth = [0.8, -0.6, 0]*util.AU_SUN
r_mars_coplanar = [-0.469, 1.45,  0]*util.AU_SUN
r_mars_incline = [-0.468, 1.44949,  00.4915]*util.AU_SUN
dt_init = (230*u.day).to(util.TU_SUN)
dt_to_burn = (115*u.day).to(util.TU_SUN)
dt_postburn = (130*u.day).to(util.TU_SUN)
meu_sun = 1*util.MEU_SUN

print("2a")
vectors =  util.get_velo_gauss_problem(r_earth, 
                                       r_mars_coplanar, 
                                       dt_init, meu_sun)["short"]
v1_init = vectors["v1"]
v2_init = vectors["v2"]

print(f"v1_init: {v1_init}")
print(f"v2_init: {v2_init}")

orbit_1 = Orbit(r_vector=r_earth, v_vector=v1_init, meu=meu_sun)
print("Initial Orbit OEs:")
print(orbit_1.orbital_elements)

print("2b")
r_at_burn, v_preburn = util.time_of_flight_universal_var(r_earth, 
                                                         v1_init, 
                                                         dt_to_burn, 
                                                         meu_sun)
print(f"Pre Burn Location: {r_at_burn}")
print(f"Pre Burn Velocity: {v_preburn.to(util.AUTU_SUN)}")

print("2c")
vectors_postburn =  util.get_velo_gauss_problem(r_at_burn,
                                                r_mars_incline, 
                                                dt_postburn, meu_sun)["short"]
v1_postburn = vectors_postburn["v1"]
v2_postburn = vectors_postburn["v2"]

print(f"v1_postburn: {v1_postburn}")
print(f"v2_postburn: {v2_postburn}")

orbit_2 = Orbit(r_vector=r_at_burn, v_vector=v1_postburn, meu=meu_sun)
print("Initial Orbit OEs:")
print(orbit_2.orbital_elements)





