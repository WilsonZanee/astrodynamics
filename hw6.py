import os

import numpy as np
from astropy import units as u

from orbits import Orbit, interplanetary_transfer_dv, OrbitalElements
import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: Earth-Jupiter interplanetary transfer ---------------------------
print("\n" + "-"*width + "\n")
print("Question 1: Earth-Jupiter interplanetary transfer\n")

mu_sun = 1*util.MEU_SUN
e_transfer = 0.7

r_earth = 1*util.AU_SUN
earth_orbit_r = 200*u.km + 1*util.DU_EARTH
mu_earth = 1*util.MEU_EARTH

r_jupiter = 5.2*util.AU_SUN
jupiter_orbit_r = 69911*u.km + 300*u.km
mu_jupiter = 1.26713e8 * (u.km**3 / u.s**2)

ra = Orbit.get_ra(r_earth, e_transfer)

# Question 1a: 
print(f"The transfer orbit with an eccentricity of 0.7 should work because\n "
      "assuming that earth is at the periapsis, the apoapsis of the oribit\n "
      f"is {ra} which is more than jupiters orbit radius of {r_jupiter}")

# Question 2a:
# Transfer Phase
transfer_energy = util.specific_energy_from_rpra(r_earth, ra, mu_sun)\
                    .to(u.km**2 / u.s**2)
vt1 = util.velo_from_energy(transfer_energy, mu_sun, r_earth)
vt2 = util.velo_from_energy(transfer_energy, mu_sun, r_jupiter)
transfer_angular_momentum = util.angular_momentum_from_periapsis(vt1, r_earth)
# Earth Escape Phase
dv1 = interplanetary_transfer_dv(earth_orbit_r,
                                 r_earth,
                                 mu_sun,
                                 transfer_angular_momentum,
                                 vt1,
                                 mu_earth, 
                                 print_v=False)
#print(dv1)
# Incoming to Jupiter
dv2 = interplanetary_transfer_dv(jupiter_orbit_r,
                                 r_jupiter,
                                 mu_sun,
                                 transfer_angular_momentum,
                                 vt2,
                                 mu_jupiter,
                                 print_v=True)
#print(dv2)
# Total dv
dv = abs(dv1) + abs(dv2)
print(dv1, dv2)
print(f"The total dv is {dv}")


# Question 2: Jupiter Fly-by --------------------------------------------------
print("\n" + "-"*width + "\n")
print("Question 2: Jupiter Fly-by\n")

escape_velo = util.get_escape_velo(mu_sun, r_jupiter)
transfer_energy = util.specific_energy_from_rpra(r_earth, r_jupiter, mu_sun)
vt1 = util.velo_from_energy(transfer_energy, mu_sun, r_earth)
vt2 = util.velo_from_energy(transfer_energy, mu_sun, r_jupiter)

# Gravity Assist
v_jupiter_mag = util.velo_from_radius(mu_sun, r_jupiter, r_jupiter)
v_jupiter = np.array([[v_jupiter_mag.value], [0], [0]])*v_jupiter_mag.unit
v_t2 = np.array([[vt2.value], [0], [0]])*vt2.unit
rp = 71398*u.km

v_final, vf_mag = util.get_gravity_assist_velo(rp, v_t2, v_jupiter, mu_jupiter, 
                                               debug=False)
print(f"To escape the Sun's influence, a S/C would need {escape_velo} of velo"
      f" where as our S/C will only have {vf_mag} of velo after getting a" 
      " gravity assist from jupiter.  Our S/C won't have enough velocity.")


# Question 3: Earth-Mars Planetary Transfer -----------------------------------
print("\n" + "-"*width + "\n")
print("Question 3: Earth-Mars Planetary Transfer\n")

r_mars = 2.28e8*u.km
mu_mars = 4.28284e4*u.km**3/u.s**2
dt = 210*u.day
r_at_earth = [-0.707, -0.707, 0]*util.AU_SUN
v_at_earth = [1.0423, -0.4124, 0]*util.AUTU_SUN
r_parked_earth = 200*u.km + 1*util.DU_EARTH
r_final_lander = 3522.2*u.km
flightpath_angle_lander = -12*u.deg
r_final_orbiter = 400*u.km + 3389.5*u.km

# 3a - Determine heliocentric r,v at Mars arrival
print("Part 3a - Spacecraft radius and velo at mars entry: ")
r_at_mars, v_at_mars = util.time_of_flight_universal_var(r_at_earth,
                                                         v_at_earth,
                                                         dt,
                                                         mu_sun)

print(f"\n {r_at_mars} = {r_at_mars.to(u.km)}"
      f"\n{v_at_mars.to(util.AUTU_SUN)} = {v_at_mars.to(u.km/u.s)}")

# 3b - Determine oe of transfer orbit
print("\nPart 3b - Transfer orbit orbital elements:")
transfer_orbit = Orbit(r_vector=r_at_earth, v_vector=v_at_earth, meu=mu_sun)
print(transfer_orbit.orbital_elements)

# 3c - Determine dv required for earth departure
print("\nPart 3c - dV for Earth departure:")
r_at_earth_scalar = np.linalg.norm(r_at_earth)
v1_transfer = np.linalg.norm(v_at_earth)
transfer_h = np.linalg.norm(np.cross(r_at_earth, v_at_earth))

dv_earth_departure = interplanetary_transfer_dv(r_parked_earth,
                                                r_at_earth_scalar,
                                                mu_sun,
                                                transfer_h,
                                                v1_transfer,
                                                mu_earth,
                                                print_v=True)
print(f"dV for Earth Departure: {dv_earth_departure}")

# 3d - determine Mars atmosphereic entry v and offset range
print(f"\nPart 3d - Mars lander and orbiter")
v2_transfer = np.linalg.norm(v_at_mars)
v_inf_mars = util.get_v_inf(mu_sun, r_mars, r_mars, transfer_h, v2_transfer)
h_lander = util.angular_momentum_enter_SOI(r_final_lander, v_inf_mars, mu_mars)
v_atmospheric_entry = (h_lander/\
                        (r_final_lander * np.cos(flightpath_angle_lander))).to(u.km/u.s)

r_mars_surf = 3380*u.km
offset_dist_surf = util.get_offset_dist(r_mars_surf, v_inf_mars, mu_mars)
offset_dist_atmo = util.get_offset_dist(r_final_lander, v_inf_mars, mu_mars)
offset_range = offset_dist_atmo - offset_dist_surf

print(f"Velo at Mars Atmospheric Entry: {v_atmospheric_entry}")
print(f"Offset range for Mars Lander: {offset_range}")

# 3e - determine dv required for mars orbiter insertion
print("\n3e - dV required for mars orbiter insertion")
dv_orbiter = interplanetary_transfer_dv(r_final_orbiter,
                                        r_mars,
                                        mu_sun,
                                        transfer_h,
                                        v2_transfer,
                                        mu_mars,
                                        print_v=True)
print(f"dV: {dv_orbiter}")
