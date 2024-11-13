import numpy as np
from numpy import cos, sin
from astropy import units as u

import two_body_util as util
import coordinate_transforms as R
from coordinate_transforms import get_local_siderail_time
from coordinate_transforms import topocentric_to_geocentrix_matrix


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
    
class RadarObservation:
    def __init__(self, range, range_rate, az, az_rate, el, el_rate, lat=None,
                 long=None, ref_zulu_theta=None, dt=None, ):
        """
        Input all angles in radians"""
        self.range = range
        self.range_rate = range_rate
        self.az = az
        self.az_rate = az_rate
        self.el = el
        self.el_rate = el_rate
        self.lat = lat
        self.long = long
        self.ref_zulu_sid_time = ref_zulu_theta
        self.dt = dt

    def get_range_vector(self):
        el = self.el.value
        az = self.az.value
        range = self.range.value
        range_vector = np.matrix([[-range*cos(el)*cos(az)],
                                  [range*cos(el)*sin(az)],
                                  [range*sin(el)]])*self.range.unit
        return range_vector
    
    def get_range_rate_vector(self):
        el = self.el.value
        el_rate = self.el_rate.value
        az = self.az.value
        az_rate = self.az_rate.value
        range = self.range.value
        range_rate = self.range_rate.value
        s_component = (-range_rate*cos(el)*cos(az) 
                        + range*sin(el)*el_rate*cos(az)
                        + range*cos(el)*sin(az)*az_rate
                        )
        e_component = (range_rate*cos(el)*sin(az) 
                        - range*sin(el)*el_rate*sin(az)
                        + range*cos(el)*cos(az)*az_rate
                        )
        z_component = (range_rate*sin(el) 
                       + range*cos(el)*el_rate 
                        )
        range_rate_vector = np.matrix([[s_component],
                                  [e_component],
                                  [z_component]])*self.range_rate.unit
        return range_rate_vector
    
    def convert_range_vector_to_radius(range_vector):
        rad_vector  = range_vector + util.earth_radius_vector 
        return rad_vector
    
    def convert_rr_IJK_to_velo(rr_IJK, radius_IJK):
        radius = np.squeeze(np.asarray(radius_IJK))*radius_IJK.unit
        omega = util.earth_rotation_velo_vector.to(u.rad/util.TU_EARTH)
        adj = np.cross(omega.value, radius.value)*rr_IJK.unit
        v = rr_IJK + adj
        return v

    def get_vectors_in_geocentric_IJK(self):
        # Get range vectors
        range_v = self.get_range_vector()
        range_rate_v = self.get_range_rate_vector()
        radius_SEZ = RadarObservation.convert_range_vector_to_radius(range_v)

        # Get D matrix
        local_sid_time = get_local_siderail_time(self.ref_zulu_sid_time, 
                                                 self.dt, 
                                                 self.long)
        D_matrix = topocentric_to_geocentrix_matrix(self.lat, local_sid_time)

        # Convert to IJK
        r = np.squeeze(np.asarray(D_matrix*radius_SEZ))*radius_SEZ.unit
        range_rate_IJK = np.squeeze(np.asarray(D_matrix*range_rate_v)
                                    )*range_rate_v.unit
        v = RadarObservation.convert_rr_IJK_to_velo(range_rate_IJK, r)

        return (r, v)

    def __str__(self):
        string = (f"Range: {self.range}\n"
                f"Range Rate: {self.range_rate}\n"
                f"Azimuth: {self.az}\n"
                f"Azimuth Rate: {self.az_rate}\n"
                f"Elevation: {self.el}\n"
                f"Elevation Rate: {self.el_rate}\n"
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
    
    def get_ra(rp, e):
        ra = ((rp.value * (1 + e)) / (1 - e))*rp.unit
        return ra

    # Finding Orbital Elements

    def get_orbital_elements_from_rv(r, v, meu):
        angular_momentum = np.cross(r,v)
        line_o_nodes = np.cross(np.array([0, 0, 1]), angular_momentum)
        e, e_vector = Orbit.get_e_from_rv(r,v, meu)

        i = Orbit.get_i(angular_momentum)
        raan = Orbit.get_raan(line_o_nodes)
        arg_periapsis = Orbit.get_arg_periapsis(line_o_nodes, e_vector)
        theta = Orbit.get_theta(r, e_vector, v)
        a = Orbit.get_a(angular_momentum, e, meu)

        oe = OrbitalElements(a.to(r.unit), e.value, i, raan, arg_periapsis, theta)
        return oe
    
    def get_a(angular_momentum, e, meu):
        if isinstance(angular_momentum, np.ndarray):
            angular_momentum = np.linalg.norm(angular_momentum)
        if isinstance(e, np.ndarray):
            e = np.linalg.norm(e)
        a = angular_momentum**2/(meu*(1-e**2))
        return a

    def get_theta(r, e_vector, v_vector):
        num = np.dot(e_vector, r)
        denom = np.linalg.norm(e_vector)*np.linalg.norm(r)
        theta = np.arccos(num/denom)
        if np.dot(r,v_vector) < 0:
            theta = 2*np.pi*u.rad - theta
        return theta

    def get_arg_periapsis(line_o_nodes, e_vector):
        num = np.dot(line_o_nodes, e_vector)
        denom = np.linalg.norm(e_vector)*np.linalg.norm(line_o_nodes)
        omega = np.arccos(num/denom)
        if e_vector[2] < 0:
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
    
def interplanetary_transfer_dv(r_parked_orbit, 
                               r_planet, 
                               mu_sun, 
                               transfer_angular_momentum,
                               v_transfer,
                               mu_planet,
                               print_v=False):
    
    v_inf = util.get_v_inf(mu_sun, 
                           r_planet, 
                           r_planet, 
                           transfer_angular_momentum,
                           v_transfer,
                           prints=print_v)
    energy_infinity = util.specific_energy_from_velo_infinity(v_inf)
    v_excess = util.velo_from_energy(energy_infinity, 
                                                mu_planet, 
                                                r_parked_orbit)
    v_planet_parked = util.velo_from_radius(mu_planet, 
                                            r_parked_orbit, 
                                            r_parked_orbit)
    if print_v:
        print(f"V excess: {v_excess.to(u.km/u.s)}")
        print(f"V parked: {v_planet_parked.to(u.km/u.s)}")


    dv = v_excess - v_planet_parked
    return dv
