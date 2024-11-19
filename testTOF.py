import numpy as np
import astropy.units as u

from orbits import Orbit, OrbitalElements, interplanetary_transfer_dv
import two_body_util as util

TOF = False
gauss = False
transfer = False
grav_assist = False
moon = True

def is_close(val_to_check, correct_val, margin=1e-3):
    vals = [val_to_check, correct_val]
    new_vals = []
    for val in vals:
        if isinstance(val, u.Quantity):
            val = val.value
        new_vals.append(val)
    if abs(new_vals[0] - new_vals[1]) < margin:
        close = True
    else:
        close = False
    return close

if TOF:
    # TOF Kepler
    print("TOF Kepler")
    e = 0.85
    a = 15000*u.km
    theta0 = 15*u.deg
    tof = 8*u.hr
    meu = 1*util.MEU_EARTH

    theta_final = util.predict_location(e, a, theta0, tof, 1, meu)
    print(theta_final)

    print(util.time_of_flight_kepler(e, a, theta0, theta_final, meu, 1).to(u.hr))

    print("TOF Kepler")
    e = 0.12
    a = 36225*u.km
    theta0 = 268*u.deg
    tof = 1.3*u.hr
    meu = 1*util.MEU_EARTH

    theta_final = util.predict_location(e, a, theta0, tof, 0, meu)
    print(theta_final)

    print(util.time_of_flight_kepler(e, a, theta0, theta_final, meu, 0).to(u.hr))

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

if gauss:
    # Gauss Problem
    # From Lecture
    r1 = [1, 0, 0]*util.DU_EARTH
    r2 = [2, 3, -1]*util.DU_EARTH
    dt = 20*util.TU_EARTH
    meu = 1*util.MEU_EARTH

    results = util.get_velo_gauss_problem(r1, r2, dt, meu, zguess=5 )
    print(results)

    correct_ans = {"short": {"v1": [1.126, 0.555, -0.185],
                            "v2": [-0.318, -0.199, 0.0664]},
                "long": {"v1": [-0.899, -0.847, 0.282],
                            "v2": [0.0477, -0.352, 0.177]}}

    passing = True
    for sol_type, vectors in results.items():
        for vector in vectors:
            for i in range(len(vector)):
                close = is_close(results[sol_type][vector][i], 
                                correct_ans[sol_type][vector][i])
                if not close:
                    print(f"{results[sol_type][vector][i]} =/=" 
                        f"{correct_ans[sol_type][vector][i]}")
                    passing = False
                    break
    if not passing:
        print(f"Test failed")

    # From Example Doc

    r1 = [.566, 1.1, -0.4]*util.DU_EARTH
    r2 = [0.33, 2.11, 0]*util.DU_EARTH
    dt = 8.55*util.TU_EARTH
    meu = 1*util.MEU_EARTH

    results = util.get_velo_gauss_problem(r1, r2, dt, meu, zguess=5 )
    print(results)

    correct_ans = {"short": {"v1": [0.263844, 0.886344, -0.12724],
                            "v2": [-0.17499, -0.47812, 0.101745]},
                "long": {"v1": [0.018046, -0.89061, -0.15975],
                            "v2": [0.195142, -0.33997, -0.25212]}}

    passing = True
    for sol_type, vectors in results.items():
        for vector in vectors:
            for i in range(len(vector)):
                close = is_close(results[sol_type][vector][i], 
                                correct_ans[sol_type][vector][i])
                if not close:
                    print(f"{results[sol_type][vector][i]} =/=" 
                        f"{correct_ans[sol_type][vector][i]}")
                    passing = False
                    break
    if not passing:
        print(f"Test failed")

    # From Example Doc

    r1 = [-1.1, 0, 0.55]*util.DU_EARTH
    r2 = [0, 0, -1]*util.DU_EARTH
    dt = 8.55*util.TU_EARTH
    meu = 1*util.MEU_EARTH

    results = util.get_velo_gauss_problem(r1, r2, dt, meu, zguess=5 )
    print(results)

    correct_ans = {"short": {"v1": [-0.9095, 0, -0.301],
                            "v2": [0.831332, 0, 0.774892]},
                "long": {"v1": [0.16025, 0, 0.943095],
                            "v2": [-1.12554, 0, 0.148431]}}

    passing = True
    for sol_type, vectors in results.items():
        for vector in vectors:
            for i in range(len(vector)):
                close = is_close(results[sol_type][vector][i], 
                                correct_ans[sol_type][vector][i])
                if not close:
                    print(f"{results[sol_type][vector][i]} =/=" 
                        f"{correct_ans[sol_type][vector][i]}")
                    passing = False
                    break
    if not passing:
        print(f"Test failed")

if transfer:
    mu_sun = 1*util.MEU_SUN
    r_earth = 1*util.AU_SUN
    r_mars = 2.28e8*u.km
    
    # Phase2
    transfer_energy = util.specific_energy_from_rpra(r_earth, r_mars, mu_sun)
    print(f"Transfer Energy: {transfer_energy}")

    v1 = util.velo_from_energy(transfer_energy, mu_sun, r_earth)
    v2 = util.velo_from_energy(transfer_energy, mu_sun, r_mars)
    print(f"v1 and v2:  {v1} {v2}")

    # Phase1
    mu_earth = 1*util.MEU_EARTH
    h = util.angular_momentum_from_periapsis(v1, r_earth)
    rp = 1*util.DU_EARTH + 300*u.km
    dv1 = interplanetary_transfer_dv(rp,
                                     r_earth,
                                     mu_sun,
                                     h,
                                     v1,
                                     mu_earth)
    print(f"Delta V 1: {dv1}")

    # Phase3
    mu_mars = 4.28284e4*u.km**3/u.s**2
    h = util.angular_momentum_from_periapsis(v1, r_mars)
    rp2 = 3780*u.km
    dv2 = interplanetary_transfer_dv(rp2,
                                     r_mars,
                                     mu_sun,
                                     h,
                                     v2,
                                     mu_mars,
                                     print_v=True)
    print(f"Delta V2: {dv2}")

    print(f"Total dV: {dv1 + dv2}")

if grav_assist:
    v_2 = np.array([[-0.50297], [0.518378], [0]])*util.AUTU_SUN
    print(f"v_init mag: {np.linalg.norm(v_2)}")
    v_mars = np.array([[-0.6116], [0.5317], [0]])*util.AUTU_SUN
    rp = 3391*u.km
    mu_mars = 4.28284e4*u.km**3/u.s**2

    print(util.get_gravity_assist_velo(rp, v_2, v_mars, mu_mars, debug=False))

if moon:
    print(f"moon_radisu at lambda=30deg: {util.get_radius_to_moon(30*u.deg)}")
    print(f"rp: {util.calc_rp_from_inf(-6.83*u.deg, 1.644*u.km/u.s, util.LUNAR_SOI, util.MU_LUNAR)}")
    