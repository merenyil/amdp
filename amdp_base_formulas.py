"""
Description
@author: Domján Ádám, extended by László Merényi
Author  = Adam Domjan (AD), László Merényi
Name    = absolute_magnetic_data_processing.py
Version = 1.0.0
Date    = 2025 october 1.

This is a library for Absolute Magnetic Data Processing (amdp) functions.

Revision History:
AD 2020 10 20: Adding null method
AD 2020 11 11: Adding trigonometric method
?? YYYY MM DD:
?? YYYY MM DD:
"""

# import sys
# sys.path.insert(1, 'P:/dev/units/python')
# sys.path.insert(1, 'P:/dev/datacq/Automata-teodolit/python_units')
# sys.path.insert(1, 'P:/dev/datacq/Automata-teodolit/measurement_post_process/python')


# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#  GLOBAL PARAMETER DEFINITIONS
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


# ----------------------------------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------------------------------

import datetime
import numpy as np
import numpy.random as RM
import numpy.linalg as LA
# import numexpr as ne                  # Faster numpy calculation

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#  SUPPORTING FUNCTIONS, CLASSES
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


##----------------------------------------------------------------------------------------------------
# Synthetic data creation - Instrument model
# ----------------------------------------------------------------------------------------------------


class synthetic_dim(object):

    def __init__(self,
                 d=0.0,
                 e=0.0,
                 o=np.array([0.0, 0.0, 0.0]),
                 s=np.array([1.0, 1.0, 1.0]),
                 azimuth_n=None,
                 azimuth_t=None,
                 mag_sigma=np.array([1.0, 1.0, 1.0]),
                 ang_sigma=np.array([5.0, 5.0])):
        """
        Creates synthetic measurements of the DIM instrument.

        param d: (rad) is the horizontal/declination collimation error.
        param e: (rad) is the verticall/inclination collimation error.
        param o: ( nT) is [1x3] the dim sensor offset. (OR noise = offset + deterministic + stochastic)
        param s: (unitless) [1x3] is the xyz component scale factor value.
        :param azimuth_n: azimuth mark angle measured in NED/world frame in radian unit (default: random, drawn fresh per instance)
        :param azimuth_t: azimuth mark angle measured in Theodolite frame in radian unit (default: random, drawn fresh per instance)
        param mag_sigma: (nT) [1x3] draw random samples from a normal (Gaussian) distribution.
        param ang_sigma: (sec) [1x2] draw random samples from a normal (Gaussian) distribution.
        """

        if azimuth_n is None: azimuth_n = RM.uniform(0., 2. * np.pi)
        if azimuth_t is None: azimuth_t = RM.uniform(0., 2. * np.pi)

        # These are the theodolite parameters
        self.delta = d
        self.epsilon = e
        self.offset = o
        self.scale = s
        self.magne_sigma = mag_sigma  # (nT)
        self.angle_sigma = ang_sigma * np.pi / 648000.  # (rad) converting second unit to radians unit

        # Noise toggles
        self.magne_normal_noise_toggle = 0
        self.angle_normal_noise_toggle = 0
        self.angle_eccentricity_toggle = 0
        self.angle_grating_error_toggle = 0  # enables graduation mark error in the theodolite circular scales
        self.angle_offset_error_toggle = 0  # gives some angle offsets in the theodolite circular scales

        # This is the true azimuth mark horizontal circle angle in the theo frame
        self.azimuth_ned = azimuth_n  # True horizontal circle angle of azimuth mark up in the ned frame
        self.azimuth_theo = azimuth_t  # True horizontal circle angle of azimuth mark up in the theo frame

        # Earth Magnetic Field:
        self.m_ned = None  # m_ned(nT) is magnetic vector field resolved in the ned frame

        # Theodolite measurements - hpa_t (rad), vpa_t (rad), msx(nT), msy(nT), msz(nT)
        # hpa_t = horizontal plane angle in theodolite frame
        # vpa_t = vertical plane angle in theodolite frame
        self.measurands = np.array([RM.uniform(0., 2. * np.pi),  # Arbitrary start for hpa_t
                                    RM.uniform(0., 2. * np.pi),  # Arbitrary start for vpa_t
                                    float("nan"),
                                    float("nan"),
                                    float("nan")])
        # self.instrument_model()

    def add_dim_noise(self,
                      mag_normal_noise=False,
                      ang_normal_noise=False,
                      ang_eccentricity=False,
                      ang_grating_error=False,
                      ang_offset_error=False):
        """
        Enables and/or disables noises

        :param mag_normal_noise: It sets the normal (Gaussian) distribution noise of the magnetic sensor.
        :param ang_normal_noise: It sets the normal (Gaussian) distribution noise of the angles.
        :param ang_eccentricity: rotation and the grating center are not coincide.
        :param ang_grating_error: error in the graduition marks of the circular scale
        :param ang_offset_error: adds offset to the theodolite angle readings
        :return:
        """

        if mag_normal_noise:
            self.magne_normal_noise_toggle = 1
        else:
            self.magne_normal_noise_toggle = 0

        if ang_normal_noise:
            self.angle_normal_noise_toggle = 1
        else:
            self.angle_normal_noise_toggle = 0

        if ang_eccentricity:
            self.angle_eccentricity_toggle = 1
        else:
            self.angle_eccentricity_toggle = 0

        if ang_grating_error:
            self.angle_grating_error_toggle = 1
        else:
            self.angle_grating_error_toggle = 0

        if ang_offset_error:
            self.angle_offset_error_toggle = 1
        else:
            self.angle_offset_error_toggle = 0

    def set_magne_stdev(self, sigma=np.array([1., 1., 1.])):
        """
        It sets the magnetic sensor standard deviation of the normal (Gaussian) distribution.

        :param sigma: [nT] default is 1 nT st. dev.
        :return:
        """
        self.magne_sigma = sigma

    def set_angle_stdev(self, sigma=np.array([5., 5.])):
        """
        It sets the angle standard deviation of the normal (Gaussian) distribution.
        Output in arc second units.

        :param sigma: [arc seconds] default is 5" arc seconds st. dev.
        :return:
        """

        self.angle_sigma = sigma * np.pi / 648000.  # convert seconds to radians

    def add_azimuth_info(self, azimuth_n=None, azimuth_t=None):
        """
        If inputs are not given thab by default this function draw samples from a uniform distribution
        for the azimuth mark angles. Samples are uniformly distributed over the half-open interval [0, 360).

        :param azimuth_n: azimuth mark angle measured in NED frame in radian unit
        :param azimuth_t: azimuth mark angle measured in Theodolite frame in radian unit
        :return:
        """

        if azimuth_n is None: azimuth_n = RM.uniform(0., 2. * np.pi)
        if azimuth_t is None: azimuth_t = RM.uniform(0., 2. * np.pi)

        self.azimuth_ned = azimuth_n
        self.azimuth_theo = azimuth_t

    def instrument_model(self):
        """
        The declination inclination magnetometer (DIM) instrument model. It calculates the xyz components of
        the DIM's sensor based on the theodolite parameters (d,e,b,c) and the world magnetic data (F,D_ned,I_ned)
        or (X_ned, Y_ned, Z_ned).

        In the real world:
        We measure:
            F: ( nT) is Nx1 vector which is contain the magnetic total field values at time of the measurments.
                Measured by a scalar magnetometer.
            hpa_ned: (rad) is Nx1 horizontal circle (plane) angle readings vector.
                hpa_ned = (hpa_theo - azimuth_theo + azimuth_ned) % (2 * np.pi)
            vpa_ned: (rad) is Nx1 vertical circle (plane) angle readings vector.
                vpa_ned = (vpa_theo - (np.pi / 2)) % (2 * np.pi)

        We do not know: D_ned, I_ned and the theodolite parameters (d, e, b, c)
            D_ned: (rad) is the declination of the time of the measurements.
            I_ned: (rad) is the inclination of the time of the measurements.
            d: (rad) is the delta angle or horizontal collimation angle
            e: (rad) is the delta angle or horizontal collimation angle
            b: ( nT) is the offset vector
            c: (unitless) is the scale vector

        :return: dim_sensor is Nx3 matrix which contains the magnetic components in the dim/theodolite frame
        """

        # scale 'rotation'
        R_s = np.diag(self.scale)

        # delta rotation
        R_d = np.array([[np.cos(self.delta), np.sin(self.delta), 0],
                        [(-1) * np.sin(self.delta), np.cos(self.delta), 0],
                        [0, 0, 1]])

        # epsilon rotation
        R_e = np.array([[np.cos(self.epsilon), 0, np.sin(self.epsilon)],
                        [0, 1, 0],
                        [(-1) * np.sin(self.epsilon), 0, np.cos(self.epsilon)]])

        # from vpa_theo to vpa_ned
        # vpa_ned = (self.measurands[1] - (np.pi / 2)) % (2 * np.pi)
        vpa_ned = self.vpa_tn(self.measurands[1])

        # vertical circle (plane) angle rotation
        R_v = np.array([[np.cos(vpa_ned), 0, np.sin(vpa_ned)],
                        [0, 1, 0],
                        [-np.sin(vpa_ned), 0, np.cos(vpa_ned)]])

        # from hpa_theo to hpa_ned
        # hpa_ned = (self.measurands[0] - self.azimuth_theo + self.azimuth_ned) % (2 * np.pi)
        hpa_ned = self.hpa_tn(self.measurands[0])

        # horizontal circle (plane) angle rotation
        R_h = np.array([[np.cos(hpa_ned), np.sin(hpa_ned), 0],
                        [(-1) * np.sin(hpa_ned), np.cos(hpa_ned), 0],
                        [0, 0, 1]])

        self.measurands[2:5] = R_s @ R_d @ R_e @ R_v @ R_h @ self.m_ned + self.offset

    def rotate(self, input1, input2, order='hpa, vpa'):
        """
        It commands the DIM into the desired angles using PID control.

        :param input1: desired horizontal/vertical circle angle resolved in theodolite frame
        :param input2: desired vertical/horizontal circle angle resolved in theodolite frame
        :param order: defines which input is which
        :return:
        """

        # TODO - mimic theodolite rotation
        # while (abs(self.measurands[0] - hpa_theo_cmd) < tol) && (abs(self.measurands[1] - vpa_theo_cmd) < tol)

        if order in ['h, v', 'hpa, vpa', 'hca, vca', 'horizontal plane angle, verical plane angle', 'Dsa, Isa']:
            hpa_theo_cmd = input1
            vpa_theo_cmd = input2
        else:
            hpa_theo_cmd = input2
            vpa_theo_cmd = input1

        self.measurands[:2] = np.array([hpa_theo_cmd, vpa_theo_cmd])

    def measure_azimuth(self, m_ned, word='mu'):
        """
        This function measures the angles of the azimuth mark in the dim/theo frame

        :param m_ned: is the NED components of the Earth Magnetic Field
        :param word: How would you like to measure azimuth mark angle?
        :return: (rad) it gives the horizontal and vertical circuit angle measurements with noise
        """

        # TODO - mimic theodolite azimuth mark measurements

        words = ['up', 'mu', 'markup', 'face left', 'circle left', 'ml']
        if word.lower() in words:
            hpa_theo = self.azimuth_theo
            vpa_theo = 1 * np.pi / 2
        else:
            hpa_theo = (self.azimuth_theo - np.pi) % (2 * np.pi)
            vpa_theo = 3 * np.pi / 2

        self.rotate(hpa_theo, vpa_theo, order='hpa, vpa')
        return self.measure(m_ned)

    def measure(self, m_ned):
        """
        It calculates the magnetic sensor components of the DIM instrument.
        It adds noise to measured angles and magnetic field to simulate real world scenario.

        :param m_ned: is the NED componets of the Earth Magnetic Field
        :return: angle and magnetic measurments corrupted with noise
        """

        self.m_ned = m_ned.flatten()
        self.instrument_model()

        angle_noise = np.zeros(2)  # vpa_t, hpa_t
        if self.angle_normal_noise_toggle:
            # Draw random samples from a normal (Gaussian) distribution for each chanel.
            vpa_noise = (RM.normal(0, self.angle_sigma[1])) % (2. * np.pi)
            hpa_noise = (RM.normal(0, self.angle_sigma[0])) % (2. * np.pi)
            angle_noise = angle_noise + np.r_[vpa_noise, hpa_noise]

        if self.angle_eccentricity_toggle:
            # Add eccentricity error to the theodolite angles
            ampli = np.array([0., 11.]) * np.pi / 648000.  # converting amplitude from arcsecond to radian unit
            phase = np.array([0., 71.]) * np.pi / 180.  # converting phase from from degree to radian unit
            angle_noise = angle_noise + ampli * np.sin(np.flip(self.measurands[0:2]) - phase)

        if self.angle_grating_error_toggle:
            # Add grating pattern errors to the theodolite angles
            ampli = np.array([0., 23.]) * np.pi / 648000.  # converting amplitude from arcsecond to radian unit
            phase = np.array([0., 277.]) * np.pi / 180.  # converting phase from from degree to radian unit
            frequ = np.array([1., 3.])  # unitless
            angle_noise = angle_noise + ampli * np.sin(frequ * (np.flip(self.measurands[0:2]) - phase))

        if self.angle_offset_error_toggle:
            # Add angle offset error to the theodolite angles
            angle_noise = angle_noise + np.array([0, 4]) * np.pi / 180.

        angles_t = (np.flip(self.measurands[0:2]) + angle_noise) % (2. * np.pi)  # vpa_t, hpa_t

        angles_n = np.array([self.vpa_tn(angles_t[0]), self.hpa_tn(angles_t[1])]) % (2. * np.pi)  # vpa_n, hpa_n

        magne_noise = np.zeros(3)
        if self.magne_normal_noise_toggle:
            # Draw random samples from a normal (Gaussian) distribution for each chanel.
            x_noise = RM.normal(0, self.magne_sigma[0], 1)  # RM.normal( mean, stdev, length)
            y_noise = RM.normal(0, self.magne_sigma[1], 1)  # RM.normal( mean, stdev, length)
            z_noise = RM.normal(0, self.magne_sigma[2], 1)  # RM.normal( mean, stdev, length)
            magne_noise = magne_noise + np.r_[x_noise, y_noise, z_noise]

        m_s = self.measurands[2:5] + magne_noise

        return angles_n, angles_t, m_s

    def angles_of_null_method(self, F_sclr, I_ned, D_ned, sequence=None, theo_angle_unit='deg'):
        """
        It calculates the angles of the null method schema for given F, I, D values

        :param F_sclr: (nT) [size=(1,1)] Total Magnetic field
        :param I_ned: (rad) [size=(1,1)] Magnetic Inclination in NED frame
        :param D_ned: (rad) [size=(1,1)] Magnetic Declination in NED frame
        :param sequence: the applied schema of null method
        :param theo_angle_unit: the synthetic 'deg', or 'gon'
        :return: (rad), size=(N,4), [c1=v_n,c2=h_n,c3=v_t,c4=h_t] the synt/true/calc angles for the null method schema
        """

        sequence = sequence if sequence else ['ED', 'WD', 'EU', 'WU', 'NU', 'SD', 'ND', 'SU', 'CH']

        def angles_of_zero(a1, a2, a3):
            dmmy = np.sqrt(a1 ** 2 + a2 ** 2 - a3 ** 2)
            s1 = np.arctan2(a2 * a3 + a1 * dmmy, a1 * a3 - a2 * dmmy) - np.pi
            s2 = - 2 * np.arctan(-a2 / a1) - s1
            return s1, s2

        vpa_n = np.empty(len(sequence))
        hpa_n = np.empty(len(sequence))

        for idx, orientation in enumerate(sequence, start=0):

            if orientation in ['EU', 'WU']:
                A1 = + self.scale[0] * F_sclr * np.cos(I_ned) * np.cos(self.delta) * np.cos(0 + self.epsilon)
                A2 = - self.scale[0] * F_sclr * np.cos(I_ned) * np.sin(self.delta)
                A3 = + self.scale[0] * F_sclr * np.sin(I_ned) * np.cos(self.delta) * np.sin(0 + self.epsilon) \
                     + self.offset[0]
                s1, s2 = angles_of_zero(A1, A2, A3)
                vpa_n[idx] = 0

                if orientation == 'WU':
                    if self.scale[0] > 0:
                        hpa_n[idx] = s1 + D_ned
                    else:
                        hpa_n[idx] = s2 + D_ned

                if orientation == 'EU':
                    if self.scale[0] > 0:
                        hpa_n[idx] = s2 + D_ned
                    else:
                        hpa_n[idx] = s1 + D_ned

            elif orientation in ['ED', 'WD']:
                A1 = + self.scale[0] * F_sclr * np.cos(I_ned) * np.cos(self.delta) * np.cos(np.pi + self.epsilon)
                A2 = - self.scale[0] * F_sclr * np.cos(I_ned) * np.sin(self.delta)
                A3 = + self.scale[0] * F_sclr * np.sin(I_ned) * np.cos(self.delta) * np.sin(np.pi + self.epsilon) \
                     + self.offset[0]
                s1, s2 = angles_of_zero(A1, A2, A3)
                vpa_n[idx] = np.pi

                if orientation == 'WD':
                    if self.scale[0] > 0:
                        hpa_n[idx] = s1 + D_ned
                    else:
                        hpa_n[idx] = s2 + D_ned

                if orientation == 'ED':
                    if self.scale[0] > 0:
                        hpa_n[idx] = s2 + D_ned
                    else:
                        hpa_n[idx] = s1 + D_ned

            elif orientation in ['NU', 'SD']:
                B1 = + self.scale[0] * F_sclr * np.cos(I_ned) * np.cos(self.delta) * (+1)
                B2 = + self.scale[0] * F_sclr * np.sin(I_ned) * np.cos(self.delta)
                B3 = - self.scale[0] * F_sclr * np.cos(I_ned) * np.sin(self.delta) * 0 + self.offset[0]
                s1, s2 = angles_of_zero(B1, B2, B3)

                if orientation == 'NU':
                    if self.scale[0] > 0:
                        vpa_n[idx] = s1 - self.epsilon
                    else:
                        vpa_n[idx] = s2 - self.epsilon

                if orientation == 'SD':
                    if self.scale[0] > 0:
                        vpa_n[idx] = s2 - self.epsilon
                    else:
                        vpa_n[idx] = s1 - self.epsilon

                hpa_n[idx] = D_ned

            elif orientation in ['ND', 'SU']:
                B1 = + self.scale[0] * F_sclr * np.cos(I_ned) * np.cos(self.delta) * (-1)
                B2 = + self.scale[0] * F_sclr * np.sin(I_ned) * np.cos(self.delta)
                B3 = - self.scale[0] * F_sclr * np.cos(I_ned) * np.sin(self.delta) * 0 + self.offset[0]
                s1, s2 = angles_of_zero(B1, B2, B3)

                if orientation == 'SU':
                    if self.scale[0] > 0:
                        vpa_n[idx] = s1 - self.epsilon
                    else:
                        vpa_n[idx] = s2 - self.epsilon

                if orientation == 'ND':
                    if self.scale[0] > 0:
                        vpa_n[idx] = s2 - self.epsilon
                    else:
                        vpa_n[idx] = s1 - self.epsilon

                hpa_n[idx] = D_ned + np.pi

            elif orientation == 'CH':

                if theo_angle_unit in ['d', 'deg', 'degree']:
                    scale_chck = 10.0 * np.pi / 10800.  # 10 min
                elif theo_angle_unit in ['g', 'gon', 'grad', 'grade', 'gradian']:
                    scale_chck = 0.2 * np.pi / 200.  # 0.2 gon
                else:
                    scale_chck = 10.0 * np.pi / 10800.  # 10 min

                vpa_n[idx] = vpa_n[idx - 1] + scale_chck if idx != 0 else I_ned  # +10 min or +0.2 gon
                hpa_n[idx] = hpa_n[idx - 1] if idx != 0 else D_ned

            elif orientation == 'MU':
                vpa_n[idx] = 0
                hpa_n[idx] = self.azimuth_ned

            elif orientation == 'MD':
                vpa_n[idx] = np.pi
                hpa_n[idx] = self.azimuth_ned + np.pi

            else:
                print('Invalid sequence entry! Please use these: [MU, MD, ED, WD, EU, WU, NU, SD, ND, SU, CH]!')
                vpa_n[idx] = np.nan
                hpa_n[idx] = np.nan

        angles_n = np.array([vpa_n, hpa_n]).T % (2. * np.pi)  # vpa_n, hpa_n
        angles_t = np.array([self.vpa_nt(vpa_n), self.hpa_nt(hpa_n)]).T % (2. * np.pi)  # vpa_t, hpa_t

        return angles_n, angles_t, sequence

    def hpa_tn(self, hpa_theo):
        """
        Converts 'Declination search angle' / 'Horizontal plane angle' from Theo to NED frame

        :param hpa_theo: (rad) 'Declination search angle' / 'Horizontal plane angle' in Theodolite frame
        :return: (rad) 'Declination search angle' in NED frame
        """
        return (hpa_theo - self.azimuth_theo + self.azimuth_ned) % (2 * np.pi)

    def hpa_nt(self, hpa_ned):
        """
        Converts 'Declination search angle' / 'Horizontal plane angle' from NED to Theo frame

        :param hpa_ned: (rad) 'Declination search angle' / 'Horizontal plane angle' in NED frame
        :return: (rad) 'Declination search angle' in Theodolite frame
        """
        return (hpa_ned - self.azimuth_ned + self.azimuth_theo) % (2 * np.pi)

    @staticmethod
    def vpa_tn(vpa_theo):
        """
        Converts 'Inclination search angle' / 'Vertical plane angle' from Theo to NED frame

        :param vpa_theo: (rad) 'Inclination search angle' / 'Horizontal plane angle' in Theodolite frame
        :return: (rad) 'Inclination search angle' in NED frame
        """

        return (vpa_theo - (np.pi / 2.)) % (2 * np.pi)

    @staticmethod
    def vpa_nt(vpa_ned):
        """
        Converts 'Inclination search angle' / 'Vertical plane angle' from NED to Theo frame

        :param vpa_ned: (rad) 'Inclination search angle' / 'Horizontal plane angle' in NED frame
        :return: (rad) 'Inclination search angle' in Theodolite frame
        """

        return (vpa_ned + (np.pi / 2.)) % (2 * np.pi)

    def get_measurands(self, theo_angle_unit='rad'):
        """
        It gives the measured values in the specific units

        :param theo_angle_unit: specify the unit of the angle output values
        :return: HPA_t, VPA_t, m_snsr_x, m_snsr_y, m_snsr_z
        """

        if theo_angle_unit in ['d', 'deg', 'degree']:
            co = 180.0 / np.pi
        elif theo_angle_unit in ['g', 'gon', 'grad', 'grade', 'gradian']:
            co = 200.0 / np.pi
        else:
            co = 1

        return np.concatenate((self.measurands[0:2] * co, self.measurands[2:]), axis=None)

    def __repr__(self, theo_angle_unit='deg'):

        if theo_angle_unit in ['d', 'deg', 'degree']:
            co = 180.0 / np.pi
        elif theo_angle_unit in ['g', 'gon', 'grad', 'grade', 'gradian']:
            co = 200.0 / np.pi
        else:
            co = 1

        line = '%8.4f, %8.4f, %+11.4f, %+11.4f, %+11.4f' % (self.measurands[0] * co, self.measurands[1] * co,
                                                            self.measurands[2], self.measurands[3], self.measurands[4])

        return line


# ----------------------------------------------------------------------------------------------------
# Function for absolute magnetic measurements - all method
# ----------------------------------------------------------------------------------------------------

def calc_doy(curr_utc_time=None):
    """
    Calculates the day of year (DoY) of a given input datetime

    param curr_utc_time: current datetime in utc unit
    return: ith day of the given year
    """

    if curr_utc_time is None: curr_utc_time = datetime.datetime.utcnow()

    begin_doy = datetime.datetime(curr_utc_time.year, 1, 1)
    begin_doy = begin_doy.toordinal()

    current_doy = datetime.datetime(curr_utc_time.year, curr_utc_time.month, curr_utc_time.day)
    current_doy = current_doy.toordinal()

    return current_doy - begin_doy + 1


def Meridian(h_ed, h_wd, h_eu, h_wu, unit='rad'):
    """
    Calculates the north and south magnetic meridian in the theodolite frame

    param h_ed: (radian) horizontal circle angle of East Down
    param h_wd: (radian) horizontal circle angle of West Down
    param h_eu: (radian) horizontal circle angle of East Up
    param h_wu: (radian) horizontal circle angle of West Up
    param unit: angle unit of input and output values, default is radian
    return: D_north_theo, D_south_theo: (radian) magnetic declination (meridian) in theo frame

    """

    if unit in ['d', 'deg', 'degree']:
        constant = 180.0
    elif unit in ['g', 'gon', 'grad', 'grade', 'gradian']:
        constant = 200.0
    else:
        constant = np.pi

    avg = np.mean([h_ed, h_wd, h_eu, h_wu])

    if (abs(h_wu - h_ed) > constant):  # Check to avoid false meridian due to epsilon/delta collimation error
        meridian = avg - constant / 2.0
    elif (abs(h_eu - h_wd) > constant):  # Check to avoid false meridian due to epsilon/delta collimation error
        meridian = avg + constant / 2.0
    elif ((h_wu - h_eu) > 0):
        meridian = avg + constant
    else:
        meridian = avg

    return meridian % (2.0 * constant), (meridian - constant) % (2.0 * constant)


# def create_angle_commands(azimuth_n_deg, azimuth_t_up_hc_tick, azimuth_t_up_vc_tick,  azimuth_t_dn_hc_tick, azimuth_t_dn_vc_tick, schema='null', latitude='N Mid Lat'):

#     import ads_instr
#     """
#     Create a schema for inverse method.

#     :param azimuth_n_deg: azimuth mark up horizontal circle angle in ned frame
#     :param azimuth_t_xx: azimuth mark up/down horizontal / vertical angles in theodolite frame
#     :param schema: Set the sequence order <(==| FOR FUTURE RELEASES
#     :param latitude: Set the latitude coordinate of the location of the measurment
#     :return: Nx2 angle commands in degree unit, 1st col = telescope angle, 2nd col = alidade angle
#                                                 1st col = vert_circ_angl, 2nd col = hori_circ_angl
#     """

#     azimuth_t_deg  = ( ads_instr.tick_to_deg(azimuth_t_up_hc_tick, False) + ads_instr.tick_to_deg(azimuth_t_dn_hc_tick, False) ) / 2.0

#     # Mark up at cícle start
#     MU = np.array( [ ads_instr.tick_to_deg(azimuth_t_up_vc_tick, True),  ads_instr.tick_to_deg(azimuth_t_up_hc_tick, False)  ] )
#     MD = np.array( [ ads_instr.tick_to_deg(azimuth_t_dn_vc_tick, True),  ads_instr.tick_to_deg(azimuth_t_dn_hc_tick, False)  ] )

#     angle_n_090 = 90.0 # (1. * np.pi / 2.)
#     angle_n_270 = 270.0 #  (3. * np.pi / 2.)
#     angle_n_360 = 360.0 # (2. * np.pi)

#     # deg2rad = np.pi / 180.

#     if   latitude == 'N High Latitude':  # North High Geographic Latitudes
#         guess = np.array([ 80., 260., 280., 100., 100. + 10. / 60.])
#     elif latitude == 'N Mid Latitude':  # North Mid Geographic Latitudes
#         guess = np.array([ 64., 244., 296., 116., 116. + 10. / 60.])
#     elif latitude == 'Equatorial':
#         guess = np.array([ 10., 190., 350., 170., 170. + 10. / 60.])
#     elif latitude == 'S Mid Latitude':
#         guess = np.array([ 62., 242., 298., 118., 118. + 10. / 60.])
#     elif latitude == 'S High Latitude':
#         guess = np.array([ 81., 261., 279.,  99.,  99. + 10. / 60.])
#     else:
#         guess = np.array([ 64., 244., 295., 116., 115. + 10. / 60.])

#     ED = np.array([angle_n_270, (angle_n_270 - azimuth_n_deg + azimuth_t_deg) % angle_n_360])  # East down
#     WD = np.array([angle_n_270, (angle_n_090 - azimuth_n_deg + azimuth_t_deg) % angle_n_360])  # West down
#     EU = np.array([angle_n_090, (angle_n_090 - azimuth_n_deg + azimuth_t_deg) % angle_n_360])  # East up
#     WU = np.array([angle_n_090, (angle_n_270 - azimuth_n_deg + azimuth_t_deg) % angle_n_360])  # West up

#     NM, SM = Meridian(ED[1], WD[1], EU[1], WU[1], 'deg')  # North and south meridian directions

#     # Iclination series
#     NU = np.array([guess[0], NM])  # North up
#     SD = np.array([guess[1], NM])  # South down
#     ND = np.array([guess[2], SM])  # North down
#     SU = np.array([guess[3], SM])  # South up
#     CH = np.array([guess[4], SM])  # (optional), Scale check (CH) measurement

#     if schema == 'null':
#         angle_commands = np.array([ MU,  MD,  ED,  WD,  EU,  WU,  NU,  SD,  ND,  SU,  CH , MU,  MD])
#         angle_commands_id =       ["MU","MD","ED","WD","EU","WU","NU","SD","ND","SU","CH","MU","MD"]
#     else:
#         angle_commands = np.array([ MU,  MD,  ED,  WD,  EU,  WU,  NU,  SD,  ND,  SU,  CH,  MU,  MD])
#         angle_commands_id =       ["MU","MD","ED","WD","EU","WU","NU","SD","ND","SU","CH","MU","MD"]

#     #for i in range(0,len(angle_commands) ):
#     #    angle_commands[i][0] = math.degrees(angle_commands[i][0])
#     #    angle_commands[i][1] = math.degrees(angle_commands[i][1])

#     return angle_commands, angle_commands_id


def hpa_tn( hpa_theo, azimuth_theo, azimuth_ned):
    """
    Converts 'Declination search angle' / 'Horizontal plane angle' from Theo to NED frame

    :param hpa_theo: (rad) 'Declination search angle' / 'Horizontal plane angle' in Theodolite frame
    :param azimuth_ned: azimuth mark angle measured in NED frame in radian unit
    :param azimuth_theo: azimuth mark angle measured in Theodolite frame in radian unit
    :return: (rad) 'Declination search angle' in NED frame
    """
    return (hpa_theo - azimuth_theo + azimuth_ned) % (2 * np.pi)


def hpa_nt( hpa_ned, azimuth_ned, azimuth_theo):
    """
    Converts 'Declination search angle' / 'Horizontal plane angle' from NED to Theo frame

    :param hpa_ned: (rad) 'Declination search angle' / 'Horizontal plane angle' in NED frame
    :param azimuth_ned: azimuth mark angle measured in NED frame in radian unit
    :param azimuth_theo: azimuth mark angle measured in Theodolite frame in radian unit
    :return: (rad) 'Declination search angle' in Theodolite frame
    """
    return (hpa_ned - azimuth_ned + azimuth_theo) % (2 * np.pi)


def vpa_tn(vpa_theo):
    """
    Converts 'Inclination search angle' / 'Vertical plane angle' from Theo to NED frame

    :param vpa_theo: (rad) 'Inclination search angle' / 'Horizontal plane angle' in Theodolite frame
    :return: (rad) 'Inclination search angle' in NED frame
    """

    return (vpa_theo - (np.pi / 2.)) % (2 * np.pi)


def vpa_nt(vpa_ned):
    """
    Converts 'Inclination search angle' / 'Vertical plane angle' from NED to Theo frame

    :param vpa_ned: (rad) 'Inclination search angle' / 'Horizontal plane angle' in NED frame
    :return: (rad) 'Inclination search angle' in Theodolite frame
    """

    return (vpa_ned + (np.pi / 2.)) % (2 * np.pi)


def statistics(data):
    """
    It gives the statistics of data column/row vector

    param data:  np.array() column or row vector
    return: [mean, std, max, min, range]
    """
    mean_val = np.mean(data)
    std_val = np.std(data)
    max_val = np.max(data)
    min_val = np.min(data)
    ptp_val = np.ptp(data)

    return np.array([mean_val, std_val, max_val, min_val, ptp_val])


# ----------------------------------------------------------------------------------------------------
# Function for absolute magnetic data processing - nonlinear least squares inverse method
# ----------------------------------------------------------------------------------------------------




def dim_instrument_model(F, vpa_n, hpa_n, I, D, e, d, o, s):
    """
    The declination inclination magnetometer instrument model.
    It calculates the x component of the sensor (sensor value) based on the model parameters
    (D, I, d, e, o, s) and the measured data (F, vca_ned, hca_ned).

    Info - Angle conversion from NED (n) to Theo (t) frame
    hpa_t = mod(hpa_n - azimuth_n + azimuth_t, 2*pi);
    We need here the zenith distance vpa_t rahter than the telescope dip angle vpa_n:[1]
    vpa_t = mod(vpa_n + pi/2, 2*pi);

    Info - Angle conversion from Theo (t) to NED (n) frame
    hpa_n = mod(hpa_t - azimuth_t + azimuth_n, 2*pi);
    We need here the telescope dip angle vpa_n rahter than the zenith distance vpa_t:[1]
    vpa_n = mod(vpa_t - pi/2, 2*pi);

    :param F: (nT) is Nx1 vector which is contain the magnetic total field
                     values at time of the measurments. Measured by a scalar magnetometer.
    :param vpa_n: (rad) is Nx1 vertical circle (plane) angle readings vector in NED frame.
    :param hpa_n: (rad) is Nx1 horizontal circle (plane) angle readings vector in NED frame.
    :param D: (rad) is the declination of the time of the measurements in NED frame.
    :param I: (rad) is the inclination of the time of the measurements in NED frame.
    :param d: (rad) is the horizontal/declination collimation error.
    :param e: (rad) is the vertical/inclination collimation error.
    :param o: (nT) is the x component of the offset (+ other noise drift, etc) of the DIM sensor
    :param s: (unitless) is the x component scale factor value
    :return: (nT) x component of the DIM sensor in the sensor frame
    """

    # For faster calculation, we pre-calculate those members that present two times in the expression:
    cos_I = np.cos(I)
    cos_delta = np.cos(d)  # 'delta', and not Declination
    vpa_n_plus_epsilon = vpa_n + e
    hpa_n_minus_decl = hpa_n - D

    a1 = + cos_I      * cos_delta     * np.cos(hpa_n_minus_decl) * np.cos(vpa_n_plus_epsilon)
    a2 = - cos_I      * np.sin(d)     * np.sin(hpa_n_minus_decl)
    a3 = + np.sin(I)  * cos_delta     * np.sin(vpa_n_plus_epsilon)

    return s * F * (a1 + a2 + a3) + o




def dim_instrument_jacobi(F, vpa_n, hpa_n, I, D, e, d, s):
    """
    The Jacobi matrix of the declination inclination magnetometer instrument model.
    This function calculates the components of the Jacobi matrix based on the model
    parameters (D,I,d,e,s) and the measured data (F, vpa_n, hpa_n).

    Info: The 'o' offset parameter has no rule here.

    :param F: ( nT) is Nx1 vector which is contain the magnetic total field values at
                    time of the measurments. Measured by a scalar magnetometer.
    :param vpa_n: (rad) is Nx1 vertical circle (plane) angle readings vector in NED/Theo frame.
    :param hpa_n: (rad) is Nx1 horizontal circle (plane) angle readings vector in NED/Theo frame.
    :param D: (rad) is the declination of the time of the measurements in NED/Theo frame.
    :param I: (rad) is the inclination of the time of the measurements in NED frame.
    :param d: (rad) is the horizontal/declination collimation error. (delta)
    :param e: (rad) is the verticall/inclination collimation error. (epsilon)
    :param s: (unitless) is the x component scale factor value
    :return: Nx6 Geometry matrix
    """

    c1 = s * F * (np.sin(vpa_n + e) * np.cos(I) * np.cos(d) -
                  np.sin(D - hpa_n) * np.sin(I) * np.sin(d) -
                  np.cos(D - hpa_n) * np.cos(vpa_n + e) * np.sin(I) * np.cos(d))
    c2 = s * F * (np.cos(D - hpa_n) * np.cos(I) * np.sin(d) -
                  np.sin(D - hpa_n) * np.cos(vpa_n + e) * np.cos(I) * np.cos(d))
    c3 = s * F * (np.cos(d) * (np.cos(vpa_n + e) * np.sin(I) - np.cos(D - hpa_n) * np.sin(vpa_n + e) * np.cos(I)))
    c4 = s * F * (np.sin(D - hpa_n) * np.cos(I) * np.cos(d) -
                  np.sin(vpa_n + e) * np.sin(I) * np.sin(d) -
                  np.cos(D - hpa_n) * np.cos(vpa_n + e) * np.cos(I) * np.sin(d))
    c5 = np.ones(np.shape(hpa_n))
    c6 = F * (np.cos(I) * np.cos(hpa_n - D) * np.cos(d) * np.cos(vpa_n + e) -
              np.cos(I) * np.sin(hpa_n - D) * np.sin(d) +
              np.sin(I) * np.sin(vpa_n + e) * np.cos(d))

    # c1 = df/dI, c2 = df/dD, c3 = df/de, c4 = df/dd, c5 = df/do, c6 = df/ds
    return np.array([c1, c2, c3, c4, c5, c6]).T


def gauss_newton(regressand: np.ndarray = None,  # Mag
                 regressors: np.ndarray = None,
                 parameters: np.ndarray = None,
                 max_iter: int = 10,
                 thresholds: np.ndarray = np.array([1e-7, 1e-3, 1e-3])
                 ):
    """
    The Gauss–Newton algorithm is used to solve non-linear least squares problems. It is a modification of Newton's
    method for finding a minimum of a function. Unlike Newton's method, the Gauss–Newton algorithm can only be used to
    minimize a sum of squared function values, but it has the advantage that second derivatives, which can be
    challenging to compute, are not required.[Wikipedia]
    Estimates the (Inclination_ned, Declination_ned, epsilon, delta, offset_x, scale_x) parameters.

    :param regressand: (msx) regressand = observed output variable of the DIM instrument model = DIM_magnetic_sensor_x
    :param regressors: (F, vpa_ned, hpa_ned) regressors = measured input vars of the DIM instrument model
    :param parameters: first for initial guess second for parameter estimations, [I, D, e, d, o, s]
    :param max_iter: Maximum allowed iteration number
    :param thresholds: thresholds for parameter, chi (RMSE), metric
    :return: parameter estimates = [I_ned (rad), D_ned (rad), e (rad), d (rad), o (nT), s (uless)]
    """

    # Initialize parameters or use guess: I, D, e, d, o_x, s_x = [0,2*pi), [0,2*pi), 0, 0, 0, +/-1
    if parameters is not None:
        parameters = parameters
    else:
        parameters = np.r_[RM.uniform(0., 2. * np.pi, 2), 0, 0, 0, np.abs(2. * RM.randint(2) - 1.)]

    # Algorithmic parameters
    stop = 0  # It controls the (while loop) iteration algorithm of the non-linear problem
    iter = 0  # iteration counter of the main while loop
    meas_nb = len(regressand)  # number of measurements
    para_nb = len(parameters)  # number of parameters
    max_iter = max_iter * para_nb  # maximum number of parameters
    dof = meas_nb - para_nb + 1  # Statistical degree of freedom

    # Print info
    message = '\nGauss-Newton algorithm'
    # print(message)

    # Main loop of the Gauss-Newton least squares fitting algorithm
    while (stop == 0) and (iter <= max_iter):

        # Increment counter
        iter += 1

        # Calculate dim magnetic sensor
        y_fit = dim_instrument_model(regressors[:, 0], regressors[:, 1], regressors[:, 2],
                                     parameters[0], parameters[1], parameters[2],
                                     parameters[3], parameters[4], parameters[5])

        # Calculate the geometry matrix
        Jacobi = dim_instrument_jacobi(regressors[:, 0], regressors[:, 1], regressors[:, 2],
                                       parameters[0], parameters[1], parameters[2],
                                       parameters[3], parameters[5])

        # Current residuals: observed minus the fitted data
        residuals = regressand - y_fit  # This is a row vector

        # chi-squared error criterion
        chi = residuals @ residuals.T

        # Auxiliary quantities
        JJ = Jacobi.T @ Jacobi

        diagJJ = np.diag(JJ)
        residualsJ = residuals @ Jacobi

        # Calculate current perturbation 'h'
        h = residualsJ @ LA.inv(JJ)

        # Update parameters
        parameters = parameters + h

        # Put angles into their ranges
        # Default: parameters[ :2] = (parameters[:2]) % (2 * np.pi)   # [0, 2*pi) range of Inclination and Declination
        parameters[0:1] = (parameters[0:1]) % (1 * np.pi) - np.pi / 2.  # [-pi/2, pi/2) range of Inclination
        parameters[1:2] = (parameters[1:2]) % (2 * np.pi)  # [0, 2*pi) range of Declination
        parameters[2:4] = (parameters[2:4] + np.pi) % (2 * np.pi) - np.pi  # [-pi, pi) range of epsilon and delta

        # Put scale factor into its range: abs(Scale) range is [0.1,10)
        # parameters[5] = np.sign(parameters[5]) * np.min([ np.max([ 0.1, np.abs(parameters[5]) ]), 10 ])

        # Convergence criteria: convergence in parameters
        if (LA.norm(np.abs(h) / (np.abs(parameters) + 10 ** (-13))) < thresholds[0]) & (iter > 2):
            # if (np.max( np.abs(h) / (np.abs(parameters) + 10**(-12)) ) < thresholds[1]) & (iter > 2):  # from [2]
            stop = 1
            info = '\nConvergence in parameters. threshold_2 = %.10f' % thresholds[0]
            # print(info)
            message += info

        # Convergence criteria: convergence in chi squared
        if (chi / dof < thresholds[1]) & (iter > 2):
            stop = 1
            info = '\nConvergence in chi squared criteria. threshold_3 = %.10f' % thresholds[1]
            # print(info)
            message += info

        # Otherwise, terminate iterations when the iteration count exceeds a pre-specified limit.
        if (iter == max_iter):
            # parameters[:] = np.NaN
            stop = 1
            info = '\nMaximum number of iteration is reached without convergence.'
            # print(info)
            message += info

    # The final residuals
    residuals = regressand - dim_instrument_model(regressors[:, 0], regressors[:, 1], regressors[:, 2],
                                                  parameters[0], parameters[1], parameters[2],
                                                  parameters[3], parameters[4], parameters[5])

    # residual statistics
    stat = statistics(residuals)

    return parameters, stat, message, residuals


# ----------------------------------------------------------------------------------------------------
# Function for absolute magnetic data processing - null method
# ----------------------------------------------------------------------------------------------------

def Declination(azimuth_ned, h_mu, h_md, h_ed, h_wd, h_eu, h_wu):
    """
    Calculates the magnetic declination in the geographic frame (North-East-Down = ned)

    param azimuth_ned: (radian) azimuth mark angle in ned frame
    param h_mu: (radian) horizontal circle angle of Mark Up position where h_mu = np.mean([h_mu_1, h_mu_2])
    param h_md: (radian) horizontal circle angle of Mark Down position where h_md = np.mean([h_md_1, h_md_2])
    param h_ed: (radian) horizontal circle angle of East Down
    param h_wd: (radian) horizontal circle angle of West Down
    param h_eu: (radian) horizontal circle angle of East Up
    param h_wu: (radian) horizontal circle angle of West Up
    :return: D_ned: (radian) magnetic declination in ned frame
    """

    # Method 1
    if (h_mu < h_md):
        h_mu_new = (h_md - np.pi) % (2.0 * np.pi)
    else:
        h_mu_new = (h_md + np.pi) % (2.0 * np.pi)

    if np.std([h_mu, h_mu_new]) > np.pi / 2.0:
        print('The zero of the horizontal circuit angle of the theodolite inside the optical collimation angle')
        azimuth_theo = float("nan")
    else:
        azimuth_theo = (h_mu + h_mu_new) / 2.0

    # Method 2
    # azimuth_theo = h_md - np.pi if (h_mu < h_md) else h_md + np.pi
    # azimuth_theo = (h_mu + azimuth_theo) / 2

    D_north_theo, _ = Meridian(h_ed, h_wd, h_eu, h_wu)

    return (azimuth_ned - azimuth_theo + D_north_theo) % (2.0 * np.pi)  # D_ned


def Inclination(v_nu, v_sd, v_nd, v_su):
    """
    Calculates the magnetic inclination in the geographic frame (North-East-Down ==> ned)

    param v_nu: (radian) vertical circle angle of North Up position
    param v_sd: (radian) vertical circle angle of South Down
    param v_nd: (radian) vertical circle angle of North Down
    param v_su: (radian) vertical circle angle of South Up
    :return: I_ned: (radian) magnetic inclination in ned frame
    """

    return (v_nu + v_sd - v_nd - v_su + 2.0 * np.pi) / 4.0  # I_ned


def Delta(h_ed, h_wd, h_eu, h_wu):
    """
    Calculates the 'delta' angle between DIM's telescope and DIM's sensor axis in the horizontal plane

    param h_ed: (radian) horizontal circle angle of East Down
    param h_wd: (radian) horizontal circle angle of West Down
    param h_eu: (radian) horizontal circle angle of East Up
    param h_wu: (radian) horizontal circle angle of West Up
    :return: delta (radian) also known as horizontal collimation error
    """

    return (h_ed + h_wd - h_eu - h_wu) / 4.0  # delta


def Epsilon_D(h_ed, h_wd, h_eu, h_wu, Inc):
    """
    Calculates the 'epsilon' angle between DIM's telescope and DIM's sensor axis in the vertical plane
    It estimates from the declination series

    param h_ed: (radian) horizontal circle angle of East Down
    param h_wd: (radian) horizontal circle angle of West Down
    param h_eu: (radian) horizontal circle angle of East Up
    param h_wu: (radian) horizontal circle angle of West Up
    param Inc: (radian) the result of Inclination(v_nu, v_sd, v_nd, v_su) function
    :return: epsilon (radian) also known as vertical collimation error
    """
    epsilon = h_eu - h_wu - h_ed + h_wd
    epsilon = epsilon - 2.0 * np.pi if (h_wd + h_eu) > (h_ed + h_wu) else epsilon + 2.0 * np.pi

    return epsilon / (4.0 * np.tan(Inc))


def Epsilon_I(v_nu, v_sd, v_nd, v_su):
    """
    Calculates the 'epsilon' angle between DIM's telescope and DIM's sensor axis in the vertical plane
    It estimates from the inclination series

    param v_nu: (radian) vertical circle angle of North Up position
    param v_sd: (radian) vertical circle angle of South Down
    param v_nd: (radian) vertical circle angle of North Down
    param v_su: (radian) vertical circle angle of South Up
    :return: epsilon (radian) also known as vertical collimation error
    """

    return np.pi - (v_nu + v_sd + v_nd + v_su) / 4.0


def Offset_D(h_ed, h_wd, h_eu, h_wu, H_D):
    """
    Calculates the offset of the DIM's sensor from the declination series

    param h_ed: (radian) horizontal circle angle of East Down
    param h_wd: (radian) horizontal circle angle of West Down
    param h_eu: (radian) horizontal circle angle of East Up
    param h_wu: (radian) horizontal circle angle of West Up
    param H_D: (nT) the average horizontal component of the magnetic field at time of the declination series where
                H_D = np.mean([F_ED, F_WD, F_EU, F_WU]) * np.cos(Inclination(v_nu, v_sd, v_nd, v_su))
    :return: offset: (nT)
    """

    return (h_eu - h_wu + h_ed - h_wd) * H_D / 4.0


def Offset_I(v_nu, v_sd, v_nd, v_su, F_I):
    """
    Calculates the offset of the DIM's sensor from the inclination series

    param v_nu: (radian) vertical circle angle of North Up position
    param v_sd: (radian) vertical circle angle of South Down
    param v_nd: (radian) vertical circle angle of North Down
    param v_su: (radian) vertical circle angle of South Up
    param F_I: (nT) is the average total magnetic field at the time of the inclination series where:
                F_I = np.mean([F_NU, F_SD, F_ND, F_SU])
    :return:
    """


    return ((v_sd + v_nd - v_su - v_nu) / 4.0 - np.pi / 2.0) * F_I


# TODO - DA - implement these functions in the future releases
# def Scale():  # need to derive its theoretical formula from [1]
# def SpotBaselineValues():  # Need scalar and variometer data




