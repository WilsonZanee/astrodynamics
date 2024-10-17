import two_body_util as util
import numpy as np
import astropy.units as u

r_init = np.array([1, 0, 1])*util.DU_EARTH
v_init = np.array([0.25, 0.5, 0])*util.DUTU_EARTH
meu = 1*util.MEU_EARTH
dt = 10*util.TU_EARTH

util.time_of_flight_universal_var(r_init, v_init, dt, meu)