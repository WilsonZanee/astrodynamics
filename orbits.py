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

        self.orbital_elements = Orbit.get_orbital_elements_from_rv(self.r, self.v)

    def get_ra(a, e):
        ra = a*(1+e)
        return ra

    def get_rp(a, e):
        ra = a*(1-e)
        return ra

    # Finding Orbital Elements

    def get_orbital_elements_from_rv(r, v, meu):
        angular_momentum = np.cross(r,v)
        line_o_nodes = np.cross(np.array([0, 0, 1]), angular_momentum)
        e, e_vector = Orbit.get_e_from_rv(r,v, meu)

        i = Orbit.get_i(angular_momentum)
        raan = Orbit.get_raan(line_o_nodes)
        arg_periapsis = Orbit.get_arg_periapsis(line_o_nodes, e_vector)
        theta = Orbit.get_theta(r, e_vector)
        a = Orbit.get_a(angular_momentum, e, meu)

        oe = OrbitalElements(a, e, i, raan, arg_periapsis, theta)
        return oe
    
    def get_a(angular_momentum, e, meu):
        if isinstance(angular_momentum, np.ndarray):
            angular_momentum = np.linalg.norm(angular_momentum)
        if isinstance(e, np.ndarray):
            e = np.linalg.norm(e)
        a = angular_momentum**2/(meu*(1-e^2))
        return a

    def get_theta(r, e_vector):
        num = np.dot(e_vector, r)
        denom = np.linalg.norm(e_vector)*np.linalg.norm(r)
        theta = np.arccos(num/denom)
        return theta

    def get_arg_periapsis(line_o_nodes, e_vector):
        num = np.dot(line_o_nodes, e_vector)
        denom = np.linalg.norm(e_vector)*np.linalg.norm(line_o_nodes)
        theta = np.arccos(num/denom)
        return theta

    def get_i(angular_momentum):
        i = np.arccos(angular_momentum[2]/np.linalg.norm(angular_momentum))
        return i

    def get_raan(line_o_nodes):
        raan = np.arccos(line_o_nodes[0]/np.linalg.norm(line_o_nodes))
        return raan

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

