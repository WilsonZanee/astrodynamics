import numpy as np
from astropy import units as u

import two_body_util as util
import coordinate_transforms as R


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

    def __str__(self):
        string = (f"Semi-Major Axis: {self.a}\n"
                f"Eccentricity: {self.e}\n"
                f"Inclination: {self.i} = {(np.rad2deg(self.i)).to(u.deg)}\n"
                f"RAAN: {self.raan} = {(np.rad2deg(self.raan)).to(u.deg)}\n"
                f"Argument of Periapsis: {self.omega} = {(np.rad2deg(self.omega)).to(u.deg)}\n"
                f"True Anomaly: {self.theta} = {(np.rad2deg(self.theta)).to(u.deg)}\n"
        )
        return string



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
        elif orbital_elements is not None:
            if self.r_vector is None and self.v_vector is None:
                self.get_attributes_oe()

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

        self.orbital_elements = Orbit.get_orbital_elements_from_rv(self.r_vector, self.v_vector, self.meu)

    def get_attributes_oe(self):
        self.p = Orbit.get_p(self.orbital_elements.e, self.orbital_elements.a) 
        self.r_mag = util.orbit_radius_from_p_eccentricity_true_anomaly(self.orbital_elements.e,
                                                                        self.p,
                                                                        self.orbital_elements.theta)
        r_perifocal = Orbit.get_r_vector_perifocal(self.r_mag, 
                                                   self.orbital_elements.theta)
        v_perifocal = Orbit.get_v_vector_perifocal(self.meu,
                                                   self.p,
                                                   self.orbital_elements.theta,
                                                   self.orbital_elements.e)
        transform_matrix = R.perifocal_to_geocentric_matrix(self.orbital_elements.raan,
                                                            self.orbital_elements.omega,
                                                            self.orbital_elements.i)
        self.r_vector = (transform_matrix*np.transpose(r_perifocal)
                                                            )*r_perifocal.unit
        self.v_vector = (transform_matrix*np.transpose(v_perifocal)
                                                            )*v_perifocal.unit


    def get_p(e, a):
        p  = a*(1-e**2)
        return p

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
        a = angular_momentum**2/(meu*(1-e**2))
        return a

    def get_theta(r, e_vector):
        num = np.dot(e_vector, r)
        denom = np.linalg.norm(e_vector)*np.linalg.norm(r)
        theta = np.arccos(num/denom)
        if num < 0:
            theta = 2*np.pi*u.rad - theta
        return theta

    def get_arg_periapsis(line_o_nodes, e_vector):
        num = np.dot(line_o_nodes, e_vector)
        denom = np.linalg.norm(e_vector)*np.linalg.norm(line_o_nodes)
        omega = np.arccos(num/denom)
        if e_vector[2] > 0:
            omega = 2*np.pi*u.rad - omega
        return omega

    def get_i(angular_momentum):
        i = np.arccos(angular_momentum[2]/np.linalg.norm(angular_momentum))
        return i

    def get_raan(line_o_nodes):
        raan = np.arccos(line_o_nodes[0]/np.linalg.norm(line_o_nodes))
        if line_o_nodes[1] < 0:
            raan = 2*np.pi*u.rad - raan
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
        e = np.linalg.norm(e_vector)

        return e, e_vector

    # Finding Radius and Velocity Vectors
    def get_r_vector_perifocal(r_mag, theta):
        r = np.matrix([r_mag.value*np.cos(theta), 
                      r_mag.value*np.sin(theta), 0])*r_mag.unit
        return r
    
    def get_v_vector_perifocal(meu, p, theta, e):
        v = np.sqrt(meu/p)*np.matrix([-np.sin(theta), (e+np.cos(theta)), 0])
        return v

