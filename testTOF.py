import numpy as np
import astropy.units as u

from orbits import Orbit, OrbitalElements
import two_body_util as util


# TOF Kepler
print("TOF Kepler")
e = 0.85
a = 15000*u.km
theta0 = 15*u.deg
tof = 8*u.hr
meu = 1*util.MEU_EARTH

theta_final = util.predict_location(e, a, theta0, tof, 1, meu)
print(theta_final)

# TOF Hyperbolic
print("TOF Hyperbolic")
r = [-0.5, 1, 0.1]*util.DU_EARTH
v = [1.2, 2, 0]*util.DUTU_EARTH
dt = 100*util.TU_EARTH
meu = 1*util.MEU_EARTH

r_final, v_final = util.time_of_flight_universal_var(r, v, dt, meu)
print(r_final, v_final)

r = [1.1, -0.4, 0]*util.DU_EARTH
v = [0.2, 0, 4.3]*util.DUTU_EARTH
dt = 54*util.TU_EARTH
meu = 1*util.MEU_EARTH

r_final, v_final = util.time_of_flight_universal_var(r, v, dt, meu)
print(r_final, v_final)

# TOF Ellipse
print("TOF Ellipse")
r = [2, 1, -1]*util.DU_EARTH
v = [0.1, 0.5, 0.2]*util.DUTU_EARTH
dt = 5*util.TU_EARTH
meu = 1*util.MEU_EARTH

r_final, v_final = util.time_of_flight_universal_var(r, v, dt, meu)
print(r_final, v_final)

r = [1, -0.22, 0]*util.DU_EARTH
v = [1, -0.1, 0.4]*util.DUTU_EARTH
dt = 25*util.TU_EARTH
meu = 1*util.MEU_EARTH

r_final, v_final = util.time_of_flight_universal_var(r, v, dt, meu)
print(r_final, v_final)