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

ra = Orbit.get_ra_fromrp(r_earth, e_transfer)

# Question 1a: 
print(f"The transfer orbit with an eccentricity of 0.7 should work because\n "
      "assuming that earth is at the periapsis, the apoapsis of the oribit\n "
      f"is {ra} which is more than jupiters orbit radius of {r_jupiter}")

# Question 1b:
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
a_mars = 1.524*util.AU_SUN
dv_orbiter = interplanetary_transfer_dv(r_final_orbiter,
                                        np.linalg.norm(r_at_mars),
                                        mu_sun,
                                        transfer_h,
                                        v2_transfer,
                                        mu_mars,
                                        a_planet=a_mars,
                                        print_v=True)
print(f"dV: {dv_orbiter}")


# Question 4: Orion S/C Lunar Orbit Insertion ---------------------------------
print("\n" + "-"*width + "\n")
print("Question 4: Orion S/C Lunar Orbit Insertion\n")

r_lunar_orbit = 252*u.km + util.LUNAR_RADIUS
v_inf = 0.5*u.km/u.s

print("4a - Epsilon_2:")
soi_lunar = util.get_SOI(
                  util.D_EARTH_LUNAR, util.LUNAR_MASS, util.MU_EARTH).to(u.km)

epsilon = util.get_epsilon2(v_inf, soi_lunar, r_lunar_orbit, util.MU_LUNAR)
print(f"Epsilon2: {epsilon.to(u.deg)}")

print("\n4b - Orbit Insertion dv and burn direction")
energy = util.specific_energy_from_velo(v_inf, util.MU_LUNAR, soi_lunar)
v_excess = util.velo_from_energy(energy, util.MU_LUNAR, r_lunar_orbit)
v_parked = util.velo_from_radius(util.MU_LUNAR, r_lunar_orbit, r_lunar_orbit)
dv = v_excess - v_parked
print(v_excess, v_parked)
print(f"dv: {dv}\nThe burn will be against the direction of motion.")


# Question 5: HLS transit to the Moon -----------------------------------------
print("\n" + "-"*width + "\n")
print("Question 5: HLS transit to the Moon\n")

print("\n5a - lambda for transfer")
r_orbit_earth = 200*u.km + 1*util.DU_EARTH
dt_transfer = 3*u.day
r_orbit_lunar = 60*u.km + util.LUNAR_RADIUS
v1 = [0.7, 0.19]*u.km/u.s
vt_lunar = np.linalg.norm(v1)
lambda1 = np.arange(42, 53, 0.01)*u.deg
best_angle = 0
best_diff = 10e6
for angle in lambda1:
      rp, epsilon2 = util.calc_rp_dv_from_patch_conic(angle, 
                                                   vt_lunar, 
                                                   r_orbit_earth)
      alt = rp - util.LUNAR_RADIUS
      diff = abs(rp - r_orbit_lunar)
      if diff.value < best_diff and epsilon2.value < 0: # e2<0 --> retrograde
            best_angle = angle
            best_diff = diff.value
            #print(f"angle: {angle}, alt: {alt}, epsilon2: {epsilon2.to(u.deg)}")
best_rp, ep2, dv = util.calc_rp_dv_from_patch_conic(best_angle, 
                                                    vt_lunar, 
                                                    r_orbit_earth,
                                                    dv=True)
best_alt = best_rp - util.LUNAR_RADIUS
print(f"HLS will achieve a lunar orbit of {best_alt} "
      f"with an lambda1 of {best_angle}")
print(f"epsilon2: {ep2.to(u.deg)}")

print("\n5a - Lunar Orbit Insertion")
print(f"dv: {dv}\nThis is a direct burn.")



      
