# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 10:23:33 2023

Software model of the ADS Theodolite.

@author: László Merényi
"""

# import datetime
import math
import numpy as np
from enum import Enum
import dim_azimuth
import angle_utils
import trigonometric_utils as tu

import misc_utils as mu

class axis_enum(Enum):
    hc = 0
    vc = 1


# _def_ADS_name = None

ENCODER_DIVISION_DICT = { "L":276960 , "K":276800, 'S':900000}  # 'S' is for simulation if a fine angle reader

M_OFFSET_OPTIM_HALF_RANGE={ "none": 0, "narrow" : 2.5,   "wide" : 145  }  # it was 40 until 2026-07-31 qqqq
M_SCALE_OPTIM_HALF_RANGE= { "none": 0, "narrow" : 0.025,  "wide" : 0.45    } 

# Global parameter defining how vc tick is converted to theo frame angle.
# It is an index of 'ADS_vctick_to_theo_frame_conv_enum'
# If "1", then vc starts at 180 theo frame, and goes until 360˘ (or 0 °) in theo frame (upside position)
# vctick_to_theo_frame_conv_type = 1

# SPEED_MONITORING_BIN_SIZE_TICKS = 100
# NOF_SPEED_MONITORING_BINS = ENCODER_DIVISION // SPEED_MONITORING_BIN_SIZE_TICKS

# speedmon_sum_tickchange = np.zeros([ 2, NOF_SPEED_MONITORING_BINS ])
# speedmon_sum_timechange = np.zeros([ 2, NOF_SPEED_MONITORING_BINS] )
# speedmon_last_bin_index = arr.array('i', [-1, -1])

def ROUND_ENCODER_TICK(f):
    return round(f)  # 'int()' or 'round()'

# def apply_encoder_type(hc_encoder_type, vc_encoder_type):
#     global HC_ENCODER_DIVISION, VC_ENCODER_DIVISION, HC_ONE_TICK_IN_SECOND, VC_ONE_TICK_IN_SECOND, hc_angle_encoder_full_range_tick_marks_pos
#     HC_ENCODER_DIVISION = ENCODER_DIVISION_DICT[hc_encoder_type]
#     VC_ENCODER_DIVISION = ENCODER_DIVISION_DICT[vc_encoder_type]
#     HC_ONE_TICK_IN_SECOND =  360.0*60.0*60.0 / HC_ENCODER_DIVISION
#     VC_ONE_TICK_IN_SECOND =  360.0*60.0*60.0 / VC_ENCODER_DIVISION
#     hc_angle_encoder_full_range_tick_marks_pos = np.arange(0, HC_ENCODER_DIVISION * 1.25, HC_ENCODER_DIVISION/4 )  # Used for plots
#     # Or [ 0, dim_instr.ENCODER_DIVISION / 4.0, dim_instr.ENCODER_DIVISION / 2.0, 3*dim_instr.ENCODER_DIVISION / 4.0,  dim_instr.ENCODER_DIVISION ]
#     vc_angle_encoder_full_range_tick_marks_pos = np.arange(0, VC_ENCODER_DIVISION * 1.25, VC_ENCODER_DIVISION/4 )  # Used for plots


# def hc_one_tick_in_rad():
#      return  pi_x_2 / ENCODER_DIVISION

# def vc_one_tick_in_rad():
#      return  pi_x_2  / ENCODER_DIVISION


# # Converts degrees into encoder ticks.
# # 'round_to_40' can minimize the encoder's error having 40 tick period.
# def deprecated_deg_to_tick(deg, deg_frame, vc_axis, restrict_0_360,  vc0 = None, round_to_40 = False):

#     if deg_frame=="raw":
#         rawtick = ROUND_ENCODER_TICK(deg / 360.0 * ENCODER_DIVISION)
#         if restrict_0_360: rawtick = rawtick % ENCODER_DIVISION
#         return rawtick
#     if deg_frame=="theo":
#         if vc_axis:
#             use_vc0(vc0)
#             deg = (deg - raw_to_theo_offset_deg - math.degrees(vc0))
#             if restrict_0_360:
#                 if deg < 0: deg += 360.0
#                 deg = deg % 360.0
#         tick = ROUND_ENCODER_TICK(deg / 360.0 * ENCODER_DIVISION)
#         if restrict_0_360: tick = (tick + ENCODER_DIVISION) % ENCODER_DIVISION
#         if (round_to_40): tick = round(tick/40) * 40
#         return tick
#     else:
#         assert False  # Not supported 'deg_frame'

# # Converts radians into encoder ticks
# def deprecated_rad_to_tick(rad, axis, restrict_0_360, vc0 = None, round_to_int = False, rad_frame = "theo" ):

#     assert isinstance(restrict_0_360, bool)
#     assert axis in ["hc", "vc"]
#     if rad_frame=="raw":
#         assert False, "Deprecated call"
#         rawtick = ROUND_ENCODER_TICK(rad / (pi_x_2) * ENCODER_DIVISION)
#         if restrict_0_360: rawtick = rawtick % ENCODER_DIVISION
#         return rawtick
#     if rad_frame=="theo":
#         if axis == "vc":
#             use_vc0(vc0)
#             rad = (rad - _vc_theo_to_tick_offset_rad - vc0)
#             if restrict_0_360:
#                 if rad < 0: rad += pi_x_2
#                 rad = rad % (pi_x_2)
#         # QQQ tick = round(rad / (pi_x_2) * ENCODER_DIVISION)
#         tick = (rad / (pi_x_2) * ENCODER_DIVISION)
#         if restrict_0_360:
#             if tick >= ENCODER_DIVISION:
#                 tick = tick - ENCODER_DIVISION
#             elif tick < 0:
#                 tick = tick + ENCODER_DIVISION
#         if round_to_int:
#             return ROUND_ENCODER_TICK(tick)
#         else:
#             return tick
#     else:
#         assert False

# def deprecated_rad_to_tick_raw(rad):
#     tick = rad / (2. * np.pi) * ENCODER_DIVISION
#     return tick


# # Converts encoder tick value into degree (Theo frame)
# def tick_to_deg(tick, vc_axis,  vc0 = None):
#     deg = (tick % ENCODER_DIVISION) * 360.0 / ENCODER_DIVISION
#     if vc_axis:
#         use_vc0(vc0)
#         deg = (deg + raw_to_theo_offset_deg + math.degrees(vc0)) % 360.0
#     return deg

# Converts encoder tick value into radians
# def tick_to_rad(tick, axis, vc0 = None):
#     assert axis in ["hc","vc"]
#     rad = tick_to_rad_raw(tick, restrict_to_0_360 = True)
#     if axis == "vc":
#         use_vc0(vc0)
#         rad = (rad + _vc_theo_to_tick_offset_rad + vc0) % pi_x_2

#     return rad


# Calculates distance between two points (two encoder ticks)
# Not yet take into account the overflow
def tick_distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))
    # return np.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))

# Takes into account the "overflow":
def tick_distance_2(a, b):
    diff = b-a
    if diff < -ENCODER_DIVISION / 2:
        diff += ENCODER_DIVISION
    elif diff > +ENCODER_DIVISION / 2:
        diff -= ENCODER_DIVISION
    return diff

# Converts a real angle (in radian) to the angle (in radian) at the encoder's resolutions
def angle_at_nearest_tick(angle, vc_axis):
    return tick_to_rad ( rad_to_tick(angle, vc_axis) )


# One single sensor with one electronics (including A/D converter, and its current configuration).
# Installation position errors on sensor holder (epsilon/delta) is also stored wtihin this class.
class mag_sensor_model_class:

    def  __init__(self, name = "M-ideal", corr_funcs = [], opt_range="none", M_scale = 1.0, M_offset = 0.0, epsilon = 0.0, delta = 0.0, M_scale_rawdata_corr = None):

        self.name = name
        self.scale  = M_scale
        # if M_offset is None: M_offset = _def_M_offset
        self.offset = M_offset

        # Position of the sensor in the sensor holder
        self.epsilon_arcmin =  angle_utils.rad_to_arcmin(epsilon)
        self.delta_arcsec =  angle_utils.rad_to_arcsec(delta)

        if name=="M-ideal":
            pass
        elif name=="ADS-MAG":
            pass
        elif name=="ADS-MAG-RevPol":
            self.scale= -1 * self.scale
        elif name=="ADS-1#2-MAG":
            self.delta_arcsec =  0  # For ADS-1#1: '-165 arcsec' before 2021.jun '109.6' after
        elif name=="SIMMAG-1":
            pass
        else:
            assert False, 'Unsupported mag sensor name: "{}"'.format(name)

        self.scale_min = self.scale - M_SCALE_OPTIM_HALF_RANGE[opt_range]
        self.scale_max = self.scale + M_SCALE_OPTIM_HALF_RANGE[opt_range]

        self.offset_min = self.offset - M_OFFSET_OPTIM_HALF_RANGE[opt_range]
        self.offset_max = self.offset + M_OFFSET_OPTIM_HALF_RANGE[opt_range]

        self.M_max = 49000         # Scale calibrated magnetometer value is allowed to be vary between -M_max ... +M_max - used as optimiztaion limits
        
        self.M_scale_rawdata_corr = M_scale_rawdata_corr

        if opt_range != 'none':

            self.epsilon_min_arcmin = -120 if opt_range == 'wide' else self.epsilon_arcmin - 1.0
            self.epsilon_max_arcmin = +120 if opt_range == 'wide' else self.epsilon_arcmin + 1.0

            self.delta_min_arcsec = -2500 if opt_range == 'wide' else self.delta_arcsec - 20
            self.delta_max_arcsec = +2500 if opt_range == 'wide' else self.delta_arcsec + 20

            # QQQ
            self.delta_min_arcsec = -8000 if opt_range == 'wide' else self.delta_arcsec - 20
            self.delta_max_arcsec = +8000 if opt_range == 'wide' else self.delta_arcsec + 20


    def epsilon(self):
        return math.radians(self.epsilon_arcmin / 60.0)

    def delta(self):
        return math.radians(self.delta_arcsec / 3600.0)

    def print(self):
        print("offset={:.2f} nT  scale={:.3f}  nT/nT  epsilon={:.1f} arcmin  delta={:.1f} arcsec".format(self.offset, self.scale, self.epsilon_arcmin, self.delta_arcsec))
        
# Used for conventional DIMs and for DS with converted angle data
class theodolite_axis_model_class:
    
    def __init__(self, parent_dim_instr, name):
        self.name = name
        self.parent_dim_instr = parent_dim_instr
        self.eccentr_error_ampl = np.nan  # If not determined, or if angle correction is not used
        self.eccentr_error_phase  = math.pi / 2.0
        self.max_eccentr_error_ampl = angle_utils.arcmin_to_rad(10)
        if name == "hc":
            self.axis = axis_enum.hc
        elif name == "vc":            
            self.axis = axis_enum.vc
        else:
            assert False
    
    
# Used for ADS with raw digitnal angle ('tick') data
class theodolite_axis_with_digital_angle_read_model_class(theodolite_axis_model_class):

    def __init__(self, parent_dim_instr, name, encoder_type):
        
        super().__init__(parent_dim_instr, name)
        self.encoder_division = ENCODER_DIVISION_DICT[encoder_type]
        self.tick_to_rad_conv = math.pi * 2 / self.encoder_division 
        self.tick_to_deg_conv = 360.0 / self.encoder_division 
        self.one_tick_in_seconds = 360.0*60.0*60.0 / self.encoder_division
        self.full_range_tick_marks_pos = np.arange(0, self.encoder_division * 1.25, self.encoder_division/4 )
        
        if self.axis == axis_enum.vc:
            self.raw_to_theo_offset_deg = 180.0 
            self.raw_to_theo_offset_rad = math.pi


    def tick_to_rad_raw(self, tick, restrict_to_0_360):
        if restrict_to_0_360:
            tick = tick % self.encoder_division
        rad = tick * self.tick_to_rad_conv
        return rad
    
    def rad_to_tick_raw(self, rad):
        tick = rad / self.tick_to_rad_conv
        return tick
    

    def tick_to_rad(self, tick, correction_functions, debug_msg):
        rad = self.tick_to_rad_raw(tick, restrict_to_0_360 = True)
        if correction_functions:
            for correction_function in correction_functions:
                if isinstance(correction_function, tu.cos_function):
                    angle_corr =  correction_function.calc(rad)
                    rad += angle_corr
                else:
                    assert correction_function.is_angle_correction()
                    rad += correction_function.calc(rad)
            rad = rad % math.tau
        
        if self.axis == axis_enum.vc:
            rad = (rad + self.raw_to_theo_offset_rad + self.parent_dim_instr.vc0) % math.tau
        return rad
    
    def deg_to_tick_raw(self, deg, restrict_0_360):
        rawtick = ROUND_ENCODER_TICK(deg / self.tick_to_deg_conv)
        if restrict_0_360: rawtick = rawtick % self.encoder_division
        return rawtick
        
    def deg_to_tick(self, deg, restrict_0_360):
        if self.axis == axis_enum.vc:
            deg = (deg - self.raw_to_theo_offset_deg - math.degrees(self.parent_dim_instr.vc0))
        return self.deg_to_tick_raw(deg, restrict_0_360)
    
    def rad_to_tick(self, rad, restrict_0_360, round_to_int = False):

        assert isinstance(restrict_0_360, bool)
        if self.axis == axis_enum.vc:
            rad = (rad - self.raw_to_theo_offset_rad - self.parent_dim_instr.vc0)
        tick =  self.rad_to_tick_raw(rad)
        if restrict_0_360:
            tick = tick % self.encoder_division
        if round_to_int:
            return ROUND_ENCODER_TICK(tick)
        else:
            return tick


# Class to model the entire ADS theodolite.
# It embeds sub-models (as classes) for one or more magnetometer sensors, and for hc
# and vc angle readers.
# As a given ADS's parameters can change over time, the period, for which the setup is created,
# must be supplied.
# 'datetime_from_to': two element list/array: 1st element is datetime-from, 2nd is datetime-to.
#  It is used to select vc0 at that time
class DIM_instrument_model_class:

    def __init__(self, config_dict,  name = None, mag_sensor_name_def = None, encoder_types = None, M_offset = 0.0, datetime_from_to = None, opt_range = "wide", vc0 = None, M_scale_rawdata_corr = None ):

        
        if config_dict is not None:
            name = config_dict.get("id", name)
            if "vc0" in config_dict.keys():
                vc0  = math.radians(config_dict.pop("vc0"))
            M_scale_rawdata_corr = config_dict.pop("M-scale-corr", M_scale_rawdata_corr)
            assert isinstance(M_scale_rawdata_corr, (float, type(None)) )
            
        assert name is not None
        if name[:3] == "ADS":
            self.type = "ADS"
        elif name[:3] == "DIM":
            self.type = "DIM"
        else:
            assert False, 'DIM instrument-id must start either with "ADS" or "DIM"'
            

        # Deprecated: convering a common angle reader setup for separate hc/vc setup
        if (encoder_types is not None) and (not isinstance(encoder_types,  dict)):
            encoder_types = {  "hc":encoder_types, "vc":encoder_types   }
            
        self.name = name # Keeping original 'name' (id), as it can be renamed below
        self.normal_sensor_updown_pos = True # Default settings for conventional DIMs

        if name == "ADS-2":
            name = "ADS-def-RSP"
            if encoder_types is None: encoder_types = { "hc": "L", "vc": "K" }
        elif name == "ADS-3":
            name = "ADS-def"
            if mag_sensor_name_def is None: mag_sensor_name_def = "ADS-MAG"
            if encoder_types is None: encoder_types = { "hc": "K", "vc": "K" }
        elif name == "ADS-4":
            name = "ADS-def-RSP"
            if mag_sensor_name_def is None: mag_sensor_name_def = "ADS-MAG"
        elif name == "ADS-5": # For the 2026 model.
            name = "ADS-def" # Note: ADS-5 virutally (not physically) had reversed sensor pos for two days in 2026: 2026-05-31..2026-06-01
        elif name in [ "ADS-VIRT","ADS-VIRT-RSP" ]:
            name = "ADS-def" if name == "ADS-VIRT" else "ADS-def-RSP"
            if mag_sensor_name_def is None: mag_sensor_name_def = "ADS-MAG"
            if encoder_types is None: encoder_types = { "hc": "S", "vc": "S" }
        elif name[:4] == "ADS-":
            name = "ADS-def-RSP" # From 2026, the default manufacturing practice is the reverse sensor pos
            
        if name in [ "ADS-def-RSP", "ADS-def" ]:
            # Generic settings for any ADS:            
            self.normal_sensor_updown_pos = (name == "ADS-def")
            if mag_sensor_name_def is None: mag_sensor_name_def = "ADS-MAG"
            if encoder_types is None:
                encoder_types = { "hc": "L", "vc": "L" } # Default is the 'Litvan'
                
        if mag_sensor_name_def is None: mag_sensor_name_def = "ADS-MAG"
                
        # Note: self.normal_sensor_updown_pos may finally be set in  'setup_for_date()'            

        self.opt_range = opt_range

        if (isinstance(mag_sensor_name_def, str)) or (not hasattr(mag_sensor_name_def, "__len__")):
            mag_sensor_name_def = [ mag_sensor_name_def ]

        M_scale_rawdata_corr = mu.ensure_list(M_scale_rawdata_corr)
        self.mag_sensors = [  ]
        for i_mag_sensor, mag_sensor_name in enumerate(mag_sensor_name_def):
            new_sensor = mag_sensor_model_class(mag_sensor_name, opt_range = opt_range, M_offset  = M_offset, M_scale_rawdata_corr = M_scale_rawdata_corr[i_mag_sensor] )
            self.mag_sensors = self.mag_sensors + [ new_sensor ]

        self.mag_measuring_range_nT = (-2000, +2000)

        self.ortho_error_arcsec = 0.0

        # if datetime_from_to == None:
        #    datetime_from_to = [ datetime.datetime(2020,1,1) , datetime.datetime.now() ]

        self.time_dependent_parameters_datetime = None
        if not datetime_from_to is None:
            self.setup_for_date(datetime_from_to)

        if opt_range != 'none':
            # self.F_pillar_diff_min = -10
            # self.F_pillar_diff_max = +10
            self.ortho_error_max_arcsec = +600 # Ortho error can vary between '-ortho_error_max_arcsec' and '+ortho_error_max_arcsec'

        match self.type:
            case "ADS":
                self.hc_axis = theodolite_axis_with_digital_angle_read_model_class(self, "hc", encoder_types["hc"])
                self.vc_axis = theodolite_axis_with_digital_angle_read_model_class(self, "vc", encoder_types["vc"])
                self.axis_angle_data_source = "digital"
                self.vc0 = vc0         # It can be None at this point: ADS raw dat files have a vc0 information in the header, that will be applied
            case "DIM":
                self.hc_axis = theodolite_axis_model_class(self, "hc")
                self.vc_axis = theodolite_axis_model_class(self, "vc")
                self.axis_angle_data_source = "manual"
                assert vc0 is None # DIM instrument does not have vc0 
                self.vc0 = 0.0
            case "_": assert False
                
        self.axis = {  axis_enum.hc : self.hc_axis, axis_enum.vc : self.vc_axis } 

        self.ext_mag_sensor_datfile_reader = None


    def hc_angle_reference_sensor_pos(self):
        return "SU" if self.normal_sensor_updown_pos else "SD"

    def setup_for_date(self,datetime_from_to):

        # for i in range(len(datetime_from_to)):
        #     if isinstance(datetime_from_to[i], datetime.datetime):
        #         datetime_from_to[i] = mu.datetime_2_exceltime(datetime_from_to[i])

        # if len(datetime_from_to)==1:
        #     datetime_from_to = [datetime_from_to[0], datetime_from_to[0]+0.1]

        # if self.name == "DIM-ideal":
        #     pass
        # elif self.name == "ADS-def":
        #     assert self.normal_sensor_updown_pos 
        # elif self.name == "ADS-def-RSP":
        #     assert  self.normal_sensor_updown_pos == False
        # elif self.name=="ADS-1#2":
        #     self.normal_sensor_updown_pos = False
        #     if datetime_from_to[0] < mu.datetime_2_exceltime(datetime.datetime(2023,3,3)):
        #         assert datetime_from_to[1] < mu.datetime_2_exceltime(datetime.datetime(2023,3,3))
        #         self.vc0 = math.radians(+2.23388887)  # Angle set by reading a horizontal MIRA.
        #     else:
        #         self.vc0 = math.radians(2.7332)   # From 2023.03.03
        # else:
        #     assert False, f'Not supported instrument: "{self.name}''

        # self.time_dependent_parameters_datetime = datetime_from_to

        # global _def_normal_sensor_updown_pos, raw_to_theo_offset_deg, _vc_theo_to_tick_offset_rad
        # _def_normal_sensor_updown_pos=self.normal_sensor_updown_pos
        # raw_to_theo_offset_deg = 180.0 if (_def_normal_sensor_updown_pos or (reverse_sensor_pos_handling_mode!=2))  else  0.0
        # assert reverse_sensor_pos_handling_mode == 4 # Added in 2025.06.04 -> only this mode should be used!!!
        # _vc_theo_to_tick_offset_rad = np.pi if (_def_normal_sensor_updown_pos or (reverse_sensor_pos_handling_mode!=2))  else  0.0
        return 

    def nof_mag_sensors(self):
        return len(self.mag_sensors)

    def get_horizontal_vc_angle(self, sensor_or_telescope_pos, unit_str):
        # assert self.time_dependent_parameters_datetime is not None
        
        if sensor_or_telescope_pos in ["FL", "FR"]:
            sensor_up_down_pos = self.FLFR_2_SUSD(sensor_or_telescope_pos)
        else:
            sensor_up_down_pos = sensor_or_telescope_pos
        
        if sensor_up_down_pos == "SU":
            if unit_str=="deg":
                return 90.0 if self.normal_sensor_updown_pos else 270.0
            elif unit_str=="rad":
               return math.pi/2.0 if self.normal_sensor_updown_pos  else 3.0*math.pi/2.0
            else:
                assert False
        elif sensor_up_down_pos == "SD":
            if unit_str=="deg":
                return 270.0 if self.normal_sensor_updown_pos  else 90.0
            elif unit_str=="rad":
               return 3.0*math.pi/2.0 if self.normal_sensor_updown_pos else math.pi/2.0
            else:
                assert False
        else:
            assert False

    # vc is in Theo frame
    # ToDo: use "get_sensor_UD_position" instead!!!!
    def is_sensor_up_at_vc(self, vc):
        # assert self.time_dependent_parameters_datetime is not None
        vc = vc % math.tau
        result = (vc < math.pi) == self.normal_sensor_updown_pos
        return result
    
    def is_telescope_in_reference_direction(self, vc_theo):
        if vc_theo % math.tau <= math.pi:
            return dim_azimuth.GEODETIC_REFERENCE_FACE_DIRECTION == "FL"
        return dim_azimuth.GEODETIC_REFERENCE_FACE_DIRECTION == "FR"
        
            
    def get_sensor_pos_at_vc(self, vc):
        return "SU" if self.is_sensor_up_at_vc( vc ) else "SD"
    
    def get_face_direction(self, vc_theo, return_format = ["FL","FR"]):  # Possible return formats: ["FL","FR"], "LR", [du.telescope_face_pos_enum.FL, du.telescope_face_pos_enum.FR]
        return return_format[0] if (vc_theo % math.tau) < math.pi else return_format[1]
        
    def get_sensor_UD_position(self, vc_theo, convert_to_default_sensor_pos):
        
        is_UP_physically = ( (vc_theo % (2*np.pi)) <= np.pi) if self.normal_sensor_updown_pos  else ( (vc_theo % (2*np.pi)) >= np.pi)

        if (not self.normal_sensor_updown_pos) and convert_to_default_sensor_pos:
            return "D" if is_UP_physically else "U"
        return "U" if is_UP_physically else "D"

    def deg_to_tick(self, deg, axis_name, restrict_0_360 = True):
        # assert  _def_normal_sensor_updown_pos == self.normal_sensor_updown_pos  # '_def_normal_sensor_updown_pos' is used by
        if axis_name == "hc":
            return self.hc.deg_to_tick(deg, restrict_0_360 )
        elif axis_name == "vc":
            return self.vc_axis.deg_to_tick(deg, restrict_0_360 )
        else:
            assert False

    def tick_to_rad(self, axis, tick, correction_functions, debug_msg):
        return self.axis[axis].tick_to_rad(tick, correction_functions, debug_msg)

    def rad_to_tick(self, axis, rad, round_to_int = False):
        if axis == axis_enum.hc:
            return self.hc.rad_to_tick(rad, round_to_int  = round_to_int )
        elif axis == axis_enum.vc:
            return self.vc_axis.rad_to_tick(rad, round_to_int  = round_to_int )
        else:
            assert False
            
    def SUSD_2_FLFR(self, SUSD):
        if SUSD == "SU"   : return "FL" if self.normal_sensor_updown_pos else "FR"
        elif SUSD == "SD" : return "FR" if self.normal_sensor_updown_pos else "FL"
        else: assert False
        
    def FLFR_2_SUSD(self, FLFR):
        if FLFR == "FL"   : return "SU" if self.normal_sensor_updown_pos else "SD"
        elif FLFR == "FR" : return "SD" if self.normal_sensor_updown_pos else "SU"
        else: assert False
        
        

# def set_def_ADS_instr(ADS_name, mag_sensor_name, encoder_type):  # ~select_ADS~
#     assert False, "Deprecated"
#     global _def_ADS_name,_def_mag_sensor_name, _def_encoder_type
#     _def_ADS_name = ADS_name
#     _def_mag_sensor_name = mag_sensor_name
#     _def_encoder_type = encoder_type
#     apply_encoder_type(encoder_type)





