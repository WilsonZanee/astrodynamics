import numpy as np

import two_body_util as util


class OrbitalElements:
    def __init__(self, semi_major_axis, eccentricity, inclination, raan, 
                 arg_periapsis, true_anomaly=None):
        self.a = semi_major_axis
        self.e = eccentricity
        self.i = inclination
        self.raan = raan
        self.omega = arg_periapsis
        self.theta = true_anomaly

    def set_true_anomaly(self, theta):
        self.theta = theta


class Orbit:
    def __init__(self, r_vector=None, v_vector=None, meu=None, 
                 orbital_elements=None):
        self.r_vector = r_vector
        self.v_vector = v_vector
        self.orbital_elements = orbital_elements
        self.meu = meu

        self.p = None
        self.angular_momentum = None
        self.spec_energy = None
        self.rp = None
        self.ra = None
        self.r_mag = None
        self.v_mag = None
        self.eccentricity = None

        if orbital_elements is None:
            if self.r_vector is not None and self.v_vector is not None:
                self.get_attributes_rv()

    def get_attributes_rv(self):
        self.angular_momentum = np.cross(self.r_vector, self.v_vector)
        self.r_mag = np.linalg.norm(self.r_vector)
        self.v_mag = np.linalg.norm(self.v_vector)
        self.spec_energy = util.specific_energy_from_velo(self.v_mag, self.meu,
                                                          self.r_mag)
        self.eccentricity = util.eccentricity_from_momentum_energy(
                                        np.linalg.norm(self.angular_momentum),
                                        self.spec_energy, self.meu)
        self.a = util.semi_major_axis_from_energy(self.spec_energy, self.meu)
        self.ra = Orbit.get_ra(self.a, self.eccentricity)
        self.rp = Orbit.get_rp(self.a, self.eccentricity)

    def get_ra(a, e):
        ra = a*(1+e)
        return ra

    def get_rp(a, e):
        ra = a*(1-e)
        return ra

    def get_oe_from_rv(r, v, meu):
        angular_momentum = np.cross(r,v)
        line_o_nodes = np.cross(np.array([0, 0, 1]), angular_momentum)

        e, e_vector = Orbit.get_e_from_rv(r,v, meu)
        a = Orbit.get_a_from_rv(r,v)
        i = Orbit.get_i_from_rv(r,v)
        raan = Orbit.get_raan_from_rv(r,v)
        omega = Orbit.get_omega_from_rv(r,v)
        theta = Orbit.get_theta_from_rv(r,v)  

        oe = OrbitalElements(a, e, i, raan, omega, theta)
        return oe

    # Finding Orbital Elements
    def get_e_from_rv(r,v, meu):
        """ Return eccentricity and eccentricity vector
        r (ndarray) - radius vector
        v (ndarray) - velocity vector
        meu (float) - Gravitational Parameter"""

        r_mag = np.linalg.norm(r)
        v_mag = np.linalg.norm(v)
        r_scaler = v_mag**2 - meu/r_mag
        v_scaler = np.dot(r, v)
        e_vector = (1/meu)*(r_scaler*r - v_scaler*v)
        e = np.linalg.norm(e)

        return e, e_vector

    def get_a_from_rv(r,v):
        return None
