import os

import numpy as np
import datetime

from astropy import units as u

from orbits import Orbit, OrbitalElements, RadarObservation
from orbits import interplanetary_transfer_dv, interplanetary_transfer_dv_vectors
import two_body_util as util

width = os.get_terminal_size().columns

# Question 1: Satellite from Radar Observation ---------------------------------
print("\n" + "-"*width + "\n")
print("Question 1: Satellite from Radar Observation\n")

print("Part a:")
mu_earth = 1*util.MEU_EARTH

range = 3.8 * util.DU_EARTH
range_rate = 0.2 * util.DUTU_EARTH
az = (18 * u.deg).to(u.rad)
az_rate = 0.01 * u.rad/util.TU_EARTH
el = (52 * u.deg).to(u.rad)
el_rate = 0.1 * u.rad/util.TU_EARTH

ref_zulu_sid_time = 20 * u.deg
lat = 45 * u.deg
long = 16 * u.deg
dt = util.get_sec("07:55:00")*u.s

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

r_sat, v_sat = obs.get_vectors_in_geocentric_IJK()

print("\nIntermediate Steps:")
print(f"{r_sat=}, {v_sat=}")
v_sat_mag = np.linalg.norm(v_sat)
r_sat_mag = np.linalg.norm(r_sat)
sat_energy = util.specific_energy_from_velo(v_sat_mag, mu_earth, r_sat_mag)
print(f"{sat_energy=}\n")

sat_orbit = Orbit(r_vector=r_sat, v_vector=v_sat, meu=mu_earth)
sat_oe_init = sat_orbit.orbital_elements
print("Satellite orbital elements at initial measurement:")
print(sat_oe_init)

print("\nPart b:")
theta_p = 360 * u.deg
sat_dt = util.time_of_flight_kepler(sat_oe_init.e, 
                                    sat_oe_init.a, 
                                    sat_oe_init.theta,
                                    theta_p,
                                    mu_earth,
                                    debug=True)
dt_pretty = datetime.timedelta(seconds=sat_dt.value)
print(f"dt to periapsis: {str(dt_pretty)} or {sat_dt:.0f}")

print("\nPart c:")

print("\nIntermediate Steps:")
r_GEO = 35786 * u.km + 1 * util.DU_EARTH

oe_periapsis = sat_oe_init
oe_periapsis.theta = theta_p

sat_orbit_peri = Orbit(orbital_elements=oe_periapsis, meu=mu_earth)
v_LEO = np.linalg.norm(sat_orbit_peri.v_vector).to(u.km/u.s)
r_LEO = np.linalg.norm(sat_orbit_peri.r_vector)

energy_LEO = util.specific_energy_from_velo(v_LEO, mu_earth, r_LEO)\
                    .to(util.DU_EARTH**2 / util.TU_EARTH**2)
v_GEO = util.velo_from_radius(mu_earth, r_GEO, r_GEO)
energy_GEO = util.specific_energy_from_velo(v_GEO, mu_earth, r_GEO)\
                    .to(util.DU_EARTH**2 / util.TU_EARTH**2)

rp_h = r_LEO
ra_h = r_GEO
energy_h = util.specific_energy_from_rpra(rp_h, ra_h, mu_earth)
vht1 = util.velo_from_energy(energy_h, mu_earth, rp_h)
vht2 = util.velo_from_energy(energy_h, mu_earth, ra_h)

# Circular LEO to hohmann elliptical
dv_sat1 = vht1 - v_LEO
# elliptical hohmann to circular GEO
dv_sat2 = dv_sat1 + (v_GEO - vht2)
print(f"{str(oe_periapsis)=}")
print(f"{energy_LEO=}")
print(f"{v_LEO=}, {v_GEO=}")
print(f"{vht1=}, {vht2=}")
print(f"{dv_sat1=}, {dv_sat2=}\n")

print(f"dv required to go to GEO: {(dv_sat1 + dv_sat2):.4f}")


# Question 2: Mars Mission ----------------------------------------------------
print("\n" + "-"*width + "\n")
print("Question 2: Mars Mission\n")
# 1 = Departure from Earth on 9/1/2037
# 2 = Arrival at Mars after 210 days
# 3 = Departure from Mars after 496 day stay
# 4 = Arrival at earth after another 210 day journey

print("Part a:")

print("\nIntermediate Steps:")
r_earth1 = [0.9422, -0.3614, -1.49e-7] * util.AU_SUN
v_earth1 = [0.3419, 0.9300, 8.543e-7] * util.AUTU_SUN
r_mars1 = [1.3871, 0.1765, -0.0304] * util.AU_SUN
v_mars1 = [-0.0714, 0.8766, 0.0201] * util.AUTU_SUN

mu_sun = 1*util.MEU_SUN

dt12 = (210 * u.day).to(util.TU_SUN)
r_mars2, v_mars2 = util.time_of_flight_universal_var(r_mars1, v_mars1, dt12, 
                                                     mu_sun)
flight_velo_vectors_12 = util.get_velo_gauss_problem(r_earth1, r_mars2, dt12,
                                                     mu_sun, 
                                                     print_convergence_table=True
                                                     )["short"]
print(f"{flight_velo_vectors_12=}")

orbit1 = Orbit(r_vector=r_earth1, v_vector=flight_velo_vectors_12["v1"],
               meu=mu_sun)
orbit2 = Orbit(r_vector=r_mars2, v_vector=flight_velo_vectors_12["v2"],
               meu=mu_sun)
print("\nOrbital Elements for Earth to Mars trajectory: ")
print(orbit1.orbital_elements)
print(f"True Anomaly at Mars: {orbit2.orbital_elements.theta.to(u.deg):.1f}")

print("Part b:")

print("\nIntermediate Steps")
dt13 = (706 * u.day).to(util.TU_SUN)
r_mars3, v_mars3 = util.time_of_flight_universal_var(r_mars1, v_mars1, dt13, 
                                                     mu_sun, debug=True)

print(f"\nMars radius at day 706: {r_mars3}")
print(f"Mars velocity at day 706: {v_mars3}")


print("\nPart c:")

print("\nIntermediate Steps")
dt34 = (210 * u.day).to(util.TU_SUN)
dt14 = (916 * u.day).to(util.TU_SUN)
r_earth4, v_earth4 = util.time_of_flight_universal_var(r_earth1, v_earth1, dt14, 
                                                     mu_sun, debug=True)

flight_velo_vectors_34 = util.get_velo_gauss_problem(r_mars3, r_earth4, dt34,
                                                     mu_sun, 
                                                     print_convergence_table=True
                                                     )["short"]
print(f"{flight_velo_vectors_34=}")

orbit3 = Orbit(r_vector=r_mars3, v_vector=flight_velo_vectors_34["v1"],
               meu=mu_sun)
orbit4 = Orbit(r_vector=r_earth4, v_vector=flight_velo_vectors_34["v2"],
               meu=mu_sun)
print("\nOrbital Elements for Mars to Earth trajectory: ")
print(orbit3.orbital_elements)
print(f"True Anomaly at Mars: {orbit4.orbital_elements.theta.to(u.deg):.1f}")

print("\nPart d:")

print("\nIntermediate Steps")
r_parked_e1 = 1000*u.km + 1*util.DU_EARTH
r_parked_m23 = 500*u.km + 3389.5*u.km
r_parked_e4 = 300*u.km + 1*util.DU_EARTH

dv1 = interplanetary_transfer_dv_vectors(r_parked_e1,
                                         flight_velo_vectors_12["v1"],
                                         v_earth1,
                                         mu_earth,
                                         print_v=True)
print(f"{dv1=}\n")

mu_mars = mu_mars = 4.28284e4*u.km**3/u.s**2
dv2 = interplanetary_transfer_dv_vectors(r_parked_m23,
                                         flight_velo_vectors_12["v2"],
                                         v_mars2,
                                         mu_mars,
                                         print_v=True)
print(f"{dv2=}\n")

dv3 = interplanetary_transfer_dv_vectors(r_parked_m23,
                                         flight_velo_vectors_34["v1"],
                                         v_mars3,
                                         mu_mars,
                                         print_v=True)
print(f"{dv3=}\n")

dv4 = interplanetary_transfer_dv_vectors(r_parked_e4,
                                         flight_velo_vectors_34["v2"],
                                         v_earth4,
                                         mu_earth,
                                         print_v=True)
print(f"{dv4=}\n")

dv_total = dv1 + dv2 + dv3 + dv4
print(f"\nTotal dv for the Mars Mission: {dv_total:.3f}")
