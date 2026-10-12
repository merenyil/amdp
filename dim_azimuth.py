# -*- coding: utf-8 -*-
"""
dim_azimuth.py
@author: Merényi László

Azimuth direction / true North calculation related functions

"""


import math
import pandas as pd
import numpy as np
import numpy.lib.recfunctions as rcf
import misc_utils as mu
import dim_utils as du
import os
import angle_utils
import geomag_utils as gu
import geomag_utils_preferences
from enum import Enum

LIMIT_FOR_FL_FR_REFMARK_READINGS_FIT__DEG = 2.0

class azimark_meas_data_source_type_enum(Enum):
    dim_dat_files      = 0   # Azimuth measured angles are together with magnetometer readings
    separate_dat_files = 1   # Azimuth angles are in separate file (ADS)
    spreadsheet_table  = 2   # Azimuth angles are in a manually edited XLS file
    
class pillar_data_source_enum(Enum):
    internal_pillar_data_table = 0
    ADS_dat_file_header        = 1
    

GEODETIC_REFERENCE_FACE_DIRECTION = "FL"  # "Face-left" / Geodetic position-I / Normal direction / (Sensor up in most DIMs)

COMPNAME__AZIMARK_INDEX = "__azi_index"

rad_angle_equality_tolerance = 1e-6

azi_meas_data_dtype = [ ('timestamp', 'f8'),  ('hc', 'f8'), ('vc', 'f8'), ('hc-tick', 'f8'), ('vc-tick', 'f8'), (COMPNAME__AZIMARK_INDEX, 'u4' )   ]    # , ('WGS84', 'f8' )

_FL_FR_angle_diff_check = None

def is_FL_FR_refmark_readings_fit(angle1, angle2):
    global _FL_FR_angle_diff_check # It can be used for error reporting, se keeping its value
    _FL_FR_angle_diff_check = abs(angle_utils.angle_distance(angle1, angle2 ) )
    _FL_FR_angle_diff_check = math.degrees(_FL_FR_angle_diff_check) 
    return abs(_FL_FR_angle_diff_check - 180.0) < LIMIT_FOR_FL_FR_REFMARK_READINGS_FIT__DEG
    

class azimark_measured_data_source_def_class:

    def __init__(self, source_type, file_name_def = None, config_dict = {}, azimark_name = None): # refmark_name = None, 
        assert isinstance(source_type, azimark_meas_data_source_type_enum)
        self.source_type = source_type
        self.file_name_def = file_name_def  # For "external-table", it is a spreadsheat file (with path), for  "separate_dat_files", it is a file name with "strftime" data masks, and with following masks : MASK_SENSORPOS: sensorpos, MASK_REFPOS: reference pos name
        assert (source_type == azimark_meas_data_source_type_enum.dim_dat_files) == (file_name_def is None) 
        self.azimark_name = azimark_name  # Required only for separate ADS DAT files, where azimuth mark name is not encoded (or encoded only in file name)
        self.max_time_diff_to_use_a_refmark_reading__days = config_dict.pop("azimuth-readings-validity-days", 1.0) # The validity lifetime of a refmark reading (in XLS file)        
        self.refmark_readings_array_index = -1 # Index of an array used in dim_dat_file_reader


# An object storing one or an associated set of azimuth direction measurement result(s).
# It typically stores azimuth measurements of one absolute measurement:
#  - in case of manual DIM meaurement: 2+2 FL+FR readings
#  - in case of ADS meaurement: 1+1 FL+FR readings
# Data can come from:
# - manual measurement XLS files
# - ADS azi readings
# - Dedicated Azi excel table (deprecated)    
# This object is then passed to 'azimuth_calculator_def_class' that manages the "theo hc angle => true angle" conversion
class complete_azimuth_measurement_data_class:

    def __init__(self):  # , static_data_source_str = '' azimark_id = "", WGS84 = None, 
         # self.azimark_id = azimark_id
         # self.WGS84 = WGS84
         # self.static_data_source_str = static_data_source_str
         # self.valid_from = self.valid_to = None
         self.angle_measurements =  { 
                             "FL": np.empty( (0),  dtype = azi_meas_data_dtype),
                             "FR": np.empty( (0),  dtype = azi_meas_data_dtype)
                        } 
         self.reset_meas_data()
      
    def reset_meas_data(self):
        self.angle_measurements["FL"].resize(0)
        self.angle_measurements["FR"].resize(0)
        assert not self.has_data("FL")
        assert not self.has_data("FR")
        self.used_azimark_indeces = set()
      
  #def reset_meas_data(self):
      # self.angle_meas_vc =  { "SU": [], "SD": []} 
      # self.angle_meas_tick_hc =  { "SU": [], "SD": []} 
      # self.angle_meas_tick_vc =  { "SU": [], "SD": []} 
      
  # def reset(self):
  #     self.refmark_name = ''
  #     self.WGS84 = None
  #     self.refpos_measured_data_source = '?'
  #     self.angle_meas_hc =  { "SU": None, "SD": None} 
  #     self.angle_meas_vc =  { "SU": None, "SD": None} 
  #     self.angle_meas_tick_hc =  { "SU": None, "SD": None} 
  #     self.angle_meas_tick_vc =  { "SU": None, "SD": None} 
  #     self.meas_timestamp = { "SU": None, "SD": None} 
  #     self.valid_from = None  # Set if data is read from external table
  #     self.valid_to   = None
      
    def use_measurement(self, FL_FR, azimark_index, timestamp, hc, vc, hc_tick = None, vc_tick = None): # , WGS84
      assert not np.isnan(hc)
      # assert not np.isnan(vc)  From 2026.10.09, we allow using nan vc angles (it is not always measrued)
      
      self.angle_measurements[FL_FR].resize( self.angle_measurements[FL_FR].shape[0] + 1 )
      new_rec = self.angle_measurements[FL_FR][-1]
      assert new_rec.base is self.angle_measurements[FL_FR]
      
      self.used_azimark_indeces.add(azimark_index)
      
      new_rec[COMPNAME__AZIMARK_INDEX] = azimark_index
      new_rec['timestamp'] = timestamp      
      # new_rec['WGS84']     = WGS84
      
      # if angle_unit == "tick":
      new_rec["hc-tick"] = hc_tick # Can be none
      new_rec["vc-tick"] = vc_tick
      new_rec["hc"] = hc
      new_rec["vc"] = vc

      # assert np.isclose(WGS84, self.WGS84, atol=rad_angle_equality_tolerance )        
      
    def get_azimark_index(self):
        if len(self.used_azimark_indeces)==1: return list(self.used_azimark_indeces)[0] 
        return None
    
      
    def print_out(self):
        #if not isinstance(self.WGS84, float):
        #if self.valid_from is not None:
        #    validity_range_str = " valid-range:{:.4f}-{:.4f}".format(self.valid_from, self.valid_to)
        #else:
        #    validity_range_str = ""
        #print("** Azimuth-meas: mark:{:} static-datasource:{:} WGS84={:.4f}° SU_hc:{:.4f}° SD_hc:{:.4f}°{}".format(self.azimark_id, self.static_data_source_str, math.degrees(self.WGS84), math.degrees(self.angle_meas_hc["SU"]), math.degrees(self.angle_meas_hc["SD"]), validity_range_str ) )
        print("Azimuth-measurement data: nr-of-FL-readings:{}  nr-of-FR-readings:{}".format(len(self.angle_measurements["FL"]), len(self.angle_measurements["FL"]) ) )
        for FLFR in ("FL","FR"):
            for i, irec in enumerate(self.angle_measurements[FLFR]):
                print(" {}. {} reading angles: date:{} hc: {:6.2f} deg  vc: {:6.2f} deg  azimark-index:{}".format(i+1, FLFR, mu.exceltime_2_datetimestr(irec["timestamp"]),  math.degrees(irec["hc"]), math.degrees(irec["vc"]), irec[COMPNAME__AZIMARK_INDEX] ) )
    
    def is_complete(self):
        return self.has_data("FL") and self.has_data("FR")
    
    def has_data(self, FL_FR):
        return len(self.angle_measurements[FL_FR]) > 0
    
    def end_of_data_feed(self, meas_data_source_str):
        self.meas_data_source_str = meas_data_source_str
        
    def first_timelabel(self):
        return min( ( np.min(self.angle_measurements["FL"]["timestamp"]), np.min(self.angle_measurements["FR"]["timestamp"]) ) )
        
    def last_timelabel(self):
        return max( ( np.max(self.angle_measurements["FL"]["timestamp"]), np.max(self.angle_measurements["FR"]["timestamp"]) ) )
    
    def nr_of_data(self, for_face_pos = '*'):
        if for_face_pos == "*": return len(self.angle_measurements["FL"]) + len(self.angle_measurements["FR"])
        return len(self.angle_measurements[for_face_pos])
    
        
  #def set_aux_params(self, meas_data_source_str, valid_from, valid_to):
  #    self.meas_data_source_str = meas_data_source_str
  #    self.valid_from = valid_from
  #    self.valid_to = valid_to
        
               


# Converts a theo hc angle into WGS84 (True) angle, using a FL/FR pair of azimuth readings 
# The pair values are typically the average readings from the start and end time of absolute measrument,
# but they can be a pair of simple readings (in case of ADS).
# It is allowed that two different (but close) azimuth marks (and azimuth directions) have been used for FL and FL measurements.
# Typical, for daytime / nithtime mark combinations in case of ADS readings at sunset/sunrise.
def h_theo_to_ned(h_theo, pos1_hc_theo, pos2_hc_theo, pos1_WGS84, pos2_WGS84, assert_result_valid = True):
    pos1_part = (pos1_hc_theo - pos1_WGS84  + math.tau) % math.tau
    pos2_part = (pos2_hc_theo - pos2_WGS84  + math.tau) % math.tau

    # FL and FR difference refers to the same WGS-direction/hc-direction combination
    distance_check = abs(angle_utils.angle_distance(pos1_part,pos2_part))
    if np.any(distance_check>0.1):
        maxindex = np.argmax(mu.ensure_list(distance_check))
        errmsg = "Error in horiz. angle calibration: SU/SD angles discrepancy: {:.2f}deg(SU) - {:.2f}deg(SD) = {:.2f}deg".format( 
            math.degrees(mu.ensure_list(pos1_part)[maxindex]), 
            math.degrees(mu.ensure_list(pos2_part)[maxindex]),
            math.degrees(max(mu.ensure_list(distance_check))) )
        assert not assert_result_valid
        return (errmsg,maxindex)

    res = h_theo - angle_utils.angle_avg(pos1_part, pos2_part)
    res = res % math.tau
    return res


# A class that manages the azimuth mark reading data (determining pillar related data).
# It has a constant data source, and then gets theodolite measurements through 'use_theo_measurement_data()'
class azimuth_calculator_def_class:

    def __init__(self, meas_position = None,  refmark_reading_related_pillar_data_source = pillar_data_source_enum.internal_pillar_data_table, references_sensor_pos = "FL", reference_telescope_position = GEODETIC_REFERENCE_FACE_DIRECTION):  # , refpos_measured_data_sources = None

        assert isinstance(refmark_reading_related_pillar_data_source, pillar_data_source_enum)
        
        assert reference_telescope_position in ["FL","FR"]
        assert references_sensor_pos in ["SU", "SD"]
        self.reference_telescope_position = reference_telescope_position
        self.reference_telescope_position_complement = "FR" if reference_telescope_position=="FL" else "FR"
        self.FLFR_2_SUSD = {  self.reference_telescope_position: references_sensor_pos, self.reference_telescope_position_complement: "SD" if references_sensor_pos == "SU" else "SU" }
        
        self.refmark_reading_related_pillar_data_source     = refmark_reading_related_pillar_data_source
        self.meas_pos = meas_position
        
        self.meas_data_source = {}  # Will be a "FL"/"FR" dictionara
        self.meas_data_timestamp = {}  # Will be a "FR"/"FLD" dictionara

        # self.refpos_const_data_source__internal_table = None  # Used if 'refmark_reading_related_pillar_data_source=="internal-table"'. Should be set in 'measurement_sites'. For format, see measurement_sites.py.
        # global _def_refpos_measured_data_source__external_table_filename
        # self.refpos_measured_data_source__external_table_filename = external_refpos_meas_table_filename # _def_refpos_measured_data_source__external_table_filename  # Used if 'refpos_measured_data_sources' has an element of "external-table". Should be set in 'measurement_sites'
        
        self.last_measured_refmark_true_angle = { "hc": { "FL": np.nan, "FR": np.nan }, "vc": { "FL": None, "FR": None }  }  # Last WGS84 azimuth belonging to last 'self.last_measured_refmark_theo_angle["hc"]["xx"]' measurement
                                                                                                                             # "vc" angle is not always know: "None" means that "vc" angle is not known.
        self.last_measured_refmark_theo_angle = { "hc": { "FL": np.nan, "FR": np.nan }, "vc": { "FL": None, "FR": None }  }  # Storing seprately the FL/FR azimuth angles measured in theo frame, to preserve the last SU/SD value for eventual missing new readings, over different day-night reference positions
                                                                                                                             # "vc" angle generally not measured ('None'), but in case of ADS, vc angle is available                       

        self.refpos_measured_data_source__external_table_pandas = None # It wil be a Panda's dataframe
        # self.last_SU_values = []  # For each 'pofpos_id', we store the last valid value
        # self.last_SD_values = []
        self.last_result = None
        # self.load_dat_file_name_list = []
        # self.load_dat = { "SU":None, "SD": None} # Data from refmark readings dat files
        
        self.new = { "FL":False, "FR": False}
        # self.act_tick_hc = { "SU":None, "SD": None}
        
        self.calculation_mode = "classical-processing-window"


    # Returns the WGS84 angle for the pier-azimark pair, Data is read out from the 'meas_pos' azimuth list.
    def get_azimark_angles(self, pillar_name, azimark_name, azi_WGS84_check = None):   # ~select_azimuth_mark~
        result = None
        if isinstance(azimark_name, bytes): azimark_name = azimark_name.decode("utf-8")
        assert self.refmark_reading_related_pillar_data_source==pillar_data_source_enum.internal_pillar_data_table
        # assert self.refpos_const_data_source__internal_table is not None, "Mesurement site is not selected, or reference mark internal table is not defined"

        # From 2025.09.22, we allow this refmark_name format: "<azimark_id>:<WGS84[deg]>
        #if ":" in azimark_id:
        #    refpos_ID, refpos_WGS84_deg = azimark_id.split(':')
        #    refpos_WGS84_deg = float(refpos_WGS84_deg)
        #    return math.radians(refpos_WGS84_deg), "data-record"
        
        
        # From 2026.08, we allow this format, if pillar is not defined:  <pillar-name>-<refmark-name>
        # act_pillar = self.meas_pos.pillars[pillar_index]
        
        hc_WGS84 = None
        
        if pillar_name == gu.PILLAR_NAME__ALL_PILLARS:  
            if "-" in azimark_name:
                act_pillar, refpos_ID = azimark_name.split("-",1)

        
        for i_record, azi_record in enumerate(self.meas_pos.azimuth_list):
            assert type(azi_record is dict)
            if (azi_record["azimark-name"]==azimark_name) and (azi_record["pier"]==pillar_name) :
                assert result is None, 'Azimuth mark ID "{}" is found more time in the table!'.format(refpos_ID)  # Ensuring that this ID exist only once
                hc_WGS84 = azi_record["hc-dir"]
                vc_angle = azi_record["vc-dir"]
                result_index = i_record
                result_info = "int-table/{}".format(result_index+1)
        
        if hc_WGS84 is  None:
            # If azi data is from Excel shert, then WGS84 can be encoded in this way: "X:10.1" Default unit is deg
            parts = azimark_name.split(":")
            assert len(parts) == 2,  'Azimuth mark ID of "{}" for pier "{}" is not found in pillars data for site \"{}\"!'.format(azimark_name, pillar_name, self.meas_pos.station_code)        
            hc_WGS84 = math.radians(float(parts[1]))
            vc_angle = None
            result_info = "from-dim-source"
        
        if azi_WGS84_check:
            assert np.isclose(hc_WGS84, azi_WGS84_check, atol = angle_utils.angle_comparison_tolerance["rad"]) 
            
        return hc_WGS84, vc_angle,  result_info

    # def get_last_value_for_refpos( self, refmark_name, SU_SD ) :
    #     if SU_SD=="SU":
    #         for i in range(len(self.last_SU_values)):
    #            if self.last_SU_values[i][0] == refmark_name:  return  self.last_SU_values[i][1]
    #         assert False, "No last value for refpos-id {}".format(refmark_name)
    #     elif SU_SD=="SD":
    #         for i in range(len(self.last_SD_values)):
    #            if self.last_SD_values[i][0] == refmark_name:  return  self.last_SD_values[i][1]
    #         assert False
    #     assert False


    def use_theo_measurement_data(self, delete_previous_measurement_data, meas_data, pillar_name, azimuth_mark_list, enable_missing_data = False):
        
        assert isinstance(meas_data, complete_azimuth_measurement_data_class)
        # assert isinstance(instr, dim_instr.DIM_instrument_model_class)
        # assert self.reference_sensor_position == instr.hc_angle_reference_sensor_pos()
        
        assert delete_previous_measurement_data, "This mode is not working"     
        assert self.calculation_mode == "classical-processing-window"
        
        azimark_index = meas_data.angle_measurements["FL"][0][COMPNAME__AZIMARK_INDEX]
        
        assert meas_data.get_azimark_index() == azimark_index # Other way to get the azimark index
        
        for FLFR in ["FL","FR"]:
            assert np.all(  meas_data.angle_measurements[FLFR][COMPNAME__AZIMARK_INDEX][:] == azimark_index), "Different asimuth marks is not handled in this version"                                                
        
        azimark_name = azimuth_mark_list[azimark_index]
        
        WGS84, azi_vc_angle, WGS84_source    = self.get_azimark_angles( pillar_name, azimark_name)
        self.nr_of_used_WGS84_dirs = 1
        self.processing_window_avg_angles = { "FL": { "hc": None, "vc":None }, "FR": { "hc": None, "vc": None } }                
        
        for FLFR in ["FL","FR"]:
            for hcvc in ["hc","vc"]: 
                # assert np.all( np.isclose(self.WGS84, meas_data.angle_measurements[FLFR][hcvc]["WGS84"][:], atol = rad_angle_equality_tolerance) )
                self.processing_window_avg_angles[FLFR][hcvc] = np.mean(  meas_data.angle_measurements[FLFR][hcvc] )                 

        if not is_FL_FR_refmark_readings_fit( self.processing_window_avg_angles["FL"]["hc"], self.processing_window_avg_angles["FR"]["hc"]):
            errmsg = "FL and FR angle mismatch: FL:{:.1f} deg  FR:{:.1f} deg  diff:{:.1f} deg".format(math.degrees(self.processing_window_avg_angles["FL"]["hc"]), math.degrees(self.processing_window_avg_angles["FR"]["hc"]), _FL_FR_angle_diff_check)
            assert False, errmsg

        self.last_measured_refmark_true_angle["hc"][self.reference_telescope_position]  = WGS84
        self.last_measured_refmark_true_angle["hc"][self.reference_telescope_position_complement]  = (WGS84 + math.pi) % math.tau
        
        self.last_measured_refmark_true_angle["vc"][self.reference_telescope_position]  = azi_vc_angle
        self.last_measured_refmark_true_angle["vc"][self.reference_telescope_position_complement]  =  None if azi_vc_angle is None else du.get_complement_vc_angle(azi_vc_angle)
        
        for FL_FR in ["FL","FR"]:
            self.new[FL_FR] = True # not np.isnan(meas_data.avg_angles[FL_FR]['hc"')
            if self.new[FL_FR]:
                 self.last_measured_refmark_theo_angle["hc"][FL_FR] = self.processing_window_avg_angles[FL_FR]["hc"]  # Storing seprately the FL/FR azimuth angles measured in theo frame, to preserve the last SU/SD value for eventual missing new readings, over different day-night reference positions
                 self.last_measured_refmark_theo_angle["vc"][FL_FR] = self.processing_window_avg_angles[FL_FR]["vc"]  # Storing seprately the FL/FR azimuth angles measured in theo frame, to preserve the last SU/SD value for eventual missing new readings, over different day-night reference positions
                 self.meas_data_source[FL_FR]    = meas_data.meas_data_source_str          # Used for debug/error check
        #         self.meas_data_timestamp[SU_SD] = meas_data.meas_timestamp[SU_SD] # Used for debug/error check
        #         self.last_measured_refmark_theo_angle["hc"][SU_SD] = meas_data.angle_meas_hc[SU_SD]
        #         if instr.hc_angle_reference_sensor_pos() == SU_SD:
        #            #  Calculating the average of the SU/SD MIRA readings
        #            self.last_measured_refmark_true_angle["hc"][SU_SD] = meas_data.WGS84
        #         else:
        #            self.last_measured_refmark_true_angle["hc"][SU_SD] = (meas_data.WGS84 + math.pi) % (math.pi*2)
        #     else:
        #         # No new SU reading: Using last value
        #         # act_MIRA_SU_hc  = self.get_last_value_for_refpos( meas_data.refmark_name, "SU" )
        #         assert enable_missing_data or (not np.isnan(self.last_measured_refmark_theo_angle["hc"][SU_SD])), "New {} reading is not valid, and has no previous SU reading. Cycle start-time:{} / enable-missing:{}".format(SU_SD, cycle_first_timelabel, enable_missing_data)

        #     self.act_tick_hc[SU_SD] = meas_data.angle_meas_tick_hc[SU_SD]
        #     if self.act_tick_hc[SU_SD] is None: self.act_tick_hc[SU_SD] = np.nan
        

        # Vertical angles are less important, just storing the newly set values (they can be NAN!)
        # self.act_SU_vc = meas_data.SU_vc
        # self.act_SD_vc = meas_data.SD_vc
        # self.act_SU_vc_tick = meas_data.SU_vc_tick
        # self.act_SD_vc_tick = meas_data.SD_vc_tick

        # Self-check: SU and SD measurements should point to the same WGS84 direction
        # if np.isfinite(self.last_measured_refmark_theo_angle["hc"]["SU"]) and np.isfinite(self.last_measured_refmark_true_angle["hc"]["SU"]) and np.isfinite(self.last_measured_refmark_theo_angle["hc"]["SD"]) and np.isfinite(self.last_measured_refmark_true_angle["hc"]["SD"]) :
        #     SU_diff = (self.last_measured_refmark_theo_angle["hc"]["SU"] - self.last_measured_refmark_true_angle["hc"]["SU"] +mu.pi_x_2) % mu.pi_x_2
        #     SD_diff = (self.last_measured_refmark_theo_angle["hc"]["SD"] - self.last_measured_refmark_true_angle["hc"]["SD"] +mu.pi_x_2) % mu.pi_x_2
        #     diff_deg = math.degrees(SU_diff-SD_diff)
        #     if diff_deg > 180: diff_deg = 360.0 - diff_deg

        #     # diff_limit_deg = 1.0
        #     # Note  kivettem, mert az XLSX feldolgozáskor élesedik az assert akkor, amikor két egymás utáni (egymástól független)
        #     # mérés követi egymást: assert abs(diff_deg) < diff_limit_deg, "MIRA-SU:{:.1f}/{:.1f}°  MIRA-SD:{:.1f}/{:.1f}° date-range:{:}-{:} diff:{:}° / limit: {:}".format(math.degrees(self.last_measured_refmark_theo_angle["hc"]["SU"]), math.degrees(self.last_measured_refmark_true_angle["hc"]["SU"]), math.degrees(self.last_measured_refmark_theo_angle["hc"]["SD"]), math.degrees(self.last_measured_refmark_true_angle["hc"]["SD"]), meas_data.valid_from, meas_data.valid_to, diff_deg, diff_limit_deg)

    # def notused__use_MIRA_data(self, start_MIRA_SU_hc, start_MIRA_SD_hc, end_MIRA_SU_hc, end_MIRA_SD_hc, pillar_azimuth_mark_ID, instr):
    #     if np.isnan(end_MIRA_SU_hc):
    #         # end MIRA data is not supplied: used when using ADS MIRA readings
    #         assert np.isnan(end_MIRA_SD_hc)
    #         self.act_MIRA_SU_hc = start_MIRA_SU_hc
    #         self.act_MIRA_SD_hc = start_MIRA_SD_hc
    #     else:
    #         self.act_MIRA_SU_hc = (start_MIRA_SU_hc+end_MIRA_SU_hc)/2.0
    #         self.act_MIRA_SD_hc = (start_MIRA_SD_hc+end_MIRA_SD_hc)/2.0

    #     if instr.hc_angle_reference_sensor_pos == "SD":
    #         temp = self.act_MIRA_SU_hc
    #         self.act_MIRA_SU_hc = self.act_MIRA_SD_hc
    #         self.act_MIRA_SD_hc = temp

    #     # Calculating the average of the SU/SD MIRA readings
    #     azimuth_theo__solution_1 =   ( self.act_MIRA_SU_hc  +  self.act_MIRA_SD_hc + ( -math.pi if (self.act_MIRA_SU_hc < self.act_MIRA_SD_hc) else  math.pi ) ) / 2.0 # From Excel

    #     # Nem biztos, hogy mindig jó implementáció!
    #     # Másik megoldás:
    #     # azimuth_theo__solution_2 =  (( (MIRA_SU_hc  + MIRA_SD_hc + np.pi) / 2.0 ) % ( math.pi ))
    #     # assert ( abs(azimuth_theo__solution_1-azimuth_theo__solution_2) < 1e-4 )

    #     self.act_azimuth_theo = azimuth_theo__solution_1
    #     self.select_azimuth_mark( pillar_azimuth_mark_ID )


    # Checks the read Excel file, if it has the required columns stored 'col_names', as a list of strings.
    def check_table_cols(self, table, source_file_name, col_names ):
        assert isinstance(table, pd.DataFrame)
        for col_name in col_names:
            if not col_name in table.columns:
                print('Column with tile "{}" is not in table of file {}!'.format(col_name, source_file_name) )
                print('Existing columns: {}'.format(table.columns))
                assert False
                
                
    #def get_act_WGS84(self):
    #    # ToDo: return a special respond, if SU and SD WGS84 origin (refmark) was different
    #    return (self.last_measured_refmark_true_angle["hc"]["SU"] + self.last_measured_refmark_true_angle["hc"]["SD"] - math.pi) % math.tau

    # def notused__use_MIRA_record(self, refpos_measured_data_source__external_table_pandas_record_index, instr, print_report_when_MIRA_changes = True):
    #       SU_col_title = 'hc SU\n(hex tick)'
    #       SD_col_title = 'hc SD\n(hex tick)'
    #     assert isinstance(self.refpos_measured_data_source__external_table_pandas, pd.DataFrame)
    #     assert refpos_measured_data_source__external_table_pandas_record_index in self.refpos_measured_data_source__external_table_pandas.index
    #     if instr.hc_angle_reference_sensor_pos == "SU":
    #     act_MIRA_SU_hc = dim_instr.tick_to_rad( int(self.refpos_measured_data_source__external_table_pandas[SU_col_title][refpos_measured_data_source__external_table_pandas_record_index], 16), False)
    #     act_MIRA_SD_hc = dim_instr.tick_to_rad( int(self.refpos_measured_data_source__external_table_pandas[SD_col_title][refpos_measured_data_source__external_table_pandas_record_index], 16), False)
    #     pillar_azimuth_mark_ID = self.refpos_measured_data_source__external_table_pandas['pillar_azimuth_mark_ID'][refpos_measured_data_source__external_table_pandas_record_index]
    #     self.use_MIRA_data(self.act_MIRA_SU_hc, self.act_MIRA_SD_hc, np.nan, np.nan, pillar_azimuth_mark_ID, instr.normal_sensor_updown_pos )
    #     if (self.last_used_MIRA_record_index != refpos_measured_data_source__external_table_pandas_record_index):
    #         self.last_used_MIRA_record_index = refpos_measured_data_source__external_table_pandas_record_index
    #         if print_report_when_MIRA_changes:
    #             print('** Using MIRA "{}" with azimuth reference "{}" **'.format(self.refpos_measured_data_source__external_table_pandas["Setup_ID"][refpos_measured_data_source__external_table_pandas_record_index], pillar_azimuth_mark_ID ))
    #             print("   MIRA-SU:{:.4f}° MIRA-SD:{:.4f}° MIRA azimuth:{:.4f}°".format(math.degrees(self.act_MIRA_SU_hc),math.degrees(self.act_MIRA_SD_hc), math.degrees(self.act_azimuth_ned)))




    # 'load_data' is data from ADS refmark data files or data from Excel files (latter not yet implemented)   
    # It must be a dict, having "FL" and "FR" keys, and a numpy array under each of these keys.
    def get_complete_azimuth_measurement_data_from_load_data(self, data_source_def, load_data, timelabel_from, timelabel_to, instr, calc_var):
        
        # Subfunction that converts ticks to rad (in case of ADS azimuth data)                
        def extract_meas_angle(FL_FR, index):
            if 'hctick' in load_data[FL_FR].dtype.names:
                assert not 'hc' in load_data[FL_FR].dtype.names
                assert instr.axis_angle_data_source == "digital"
                hc = instr.hc_axis.tick_to_rad( load_data[FL_FR][index]["hctick"], calc_var.hc_precorr_funcs, "azi-hc" )
                vc = instr.vc_axis.tick_to_rad( load_data[FL_FR][index]["vctick"], calc_var.vc_precorr_funcs, "azi-vc" )
            else:
                assert False, "Processing hc/vc angles not yet implemented"
            return hc,vc
        
        pillar_index = load_data["FL"]["pillar_index"][0]
        assert np.all(load_data["FL"]["pillar_index"] == pillar_index) 
        assert np.all(load_data["FR"]["pillar_index"] == pillar_index) 
        # WGS84, refpos_static_data_source = self.get_azimark_WGS84_angle( pillar_index, data_source_def.refmark_name)
        result = complete_azimuth_measurement_data_class() # data_source_def.refmark_name, WGS84, refpos_static_data_source
        # result.refpos_measured_data_source = "dat-file/"+data_source_def.file_name_def
        
        

        # Method until 2026.09.29 - using one FL and one FR measurement, by time checking before/after
        # FL_FR_nearest_meas = {}  # Determininf the closest azimuth measurement for the special cases
        #                        # when ADS is reinstalled on the pillar with a different orientation,
        #                        # and measurements have been performed just before and after
        #                        # This situation needs special care for pairing
        # for FL_FR in ["FL","FR"]:
        #     # Other - more simple and fast - time pairing method would be using the average (mid) timelabel
        #     nearest_from_start_index    =  (np.abs(load_data[FL_FR]['datetime'] - timelabel_from)).argmin()
        #     nearest_from_end_index      =  (np.abs(load_data[FL_FR]['datetime'] - timelabel_to)).argmin()
        #     nearest_from_start_timediff =  np.abs(load_data[FL_FR]['datetime'][nearest_from_start_index] - timelabel_from)
        #     nearest_from_end_timediff   =  np.abs(load_data[FL_FR]['datetime'][nearest_from_end_index] - timelabel_to)

        #     if nearest_from_start_timediff < nearest_from_end_timediff:
        #         time_diff = nearest_from_start_timediff
        #         nearest_record_index = nearest_from_start_index
        #     else:
        #         time_diff = nearest_from_end_timediff
        #         nearest_record_index = nearest_from_end_index
        #     if time_diff > data_source_def.max_time_diff_to_use_a_refmark_reading__days: return None,None # No SU or SD reading for this cycle
   
        #     if (not FL_FR_nearest_meas) or ( time_diff < FL_FR_nearest_meas["time_diff"]  ):
        #         FL_FR_nearest_meas["time_diff"] = time_diff
        #         FL_FR_nearest_meas["face-dir"] = FL_FR
        #         FL_FR_nearest_meas["index"] = nearest_record_index

        #     assert not np.isnan( load_data[FL_FR][nearest_record_index]["hctick"] ), "{}".format(load_data[FL_FR][nearest_record_index]["datetime"])
        #     azimark_index = 0
        #     result.use_measurement(FL_FR, azimark_index, load_data[FL_FR][nearest_record_index]['datetime'], load_data[FL_FR][nearest_record_index]["hctick"], load_data[FL_FR][nearest_record_index]["vctick"], "tick", instr )
        
        
        mid_timestamp = (timelabel_from + timelabel_to) / 2.0
        closest_measurements = {}
        for FL_FR in ["FL","FR"]:
            meas_timestamps = load_data[FL_FR]['datetime']
            is_sorted = np.all(np.diff(meas_timestamps)>= 0)  
            assert is_sorted
            closest_index = np.searchsorted(meas_timestamps, mid_timestamp)
            if (closest_index > 0):
                if closest_index == len(meas_timestamps):  # mid_timestamp is after the last azi-meas timestamp
                    closest_index -=1
                else:
                    if meas_timestamps[closest_index] - mid_timestamp > mid_timestamp - meas_timestamps[closest_index-1]:
                        closest_index -= 1

            closest_diff = abs(meas_timestamps[closest_index] - mid_timestamp) 
            if meas_timestamps[closest_index] < mid_timestamp:
                closest_side_index = "before-index"
                closest_index_before = closest_index
                if closest_index == len(meas_timestamps) - 1:
                    closest_index_after = None
                else:
                    closest_index_after = closest_index_before + 1
                    assert meas_timestamps[closest_index_after] > mid_timestamp
            else:
                closest_side_index = "after-index"
                closest_index_after  = closest_index
                if closest_index == 0:
                    closest_index_before = None
                else:
                    closest_index_before = closest_index - 1
                    assert meas_timestamps[closest_index_before] < mid_timestamp
            
            if closest_index_before is None:
                before_dir_angle = None
            else:
                before_dir_angle = instr.hc_axis.tick_to_rad( load_data[FL_FR][closest_index_before]["hctick"],  calc_var.hc_precorr_funcs, "azi-pairing" ) % math.tau
            if closest_index_after is None:
                after_dir_angle = None
            else:
                after_dir_angle = instr.hc_axis.tick_to_rad( load_data[FL_FR][closest_index_after]["hctick"],  calc_var.hc_precorr_funcs, "azi-pairing" ) % math.tau
            #if FL_FR == "FR":  
            #    if before_reference_dir_angle: before_reference_dir_angle = (before_reference_dir_angle + np.pi) % np.pi
            #    if after_reference_dir_angle:  after_reference_dir_angle  = (after_reference_dir_angle + np.pi) % np.pi
            closest_measurements[FL_FR] = { "diff": closest_diff, "closest-side-key":closest_side_index, "before-index":closest_index_before, "after-index": closest_index_after, "before-ref-angle": before_dir_angle, "after-ref-angle":after_dir_angle, "closest-ref-angle": before_dir_angle if closest_side_index == "before-index" else after_dir_angle }

               
        # We must handle the special case when ADS/DIM is reinstalled on the pillar with a different orientation,
        # and non conventional ordered measurements have been performed just before and after.
        # This situation needs special care for pairing!
        # First, checking if closest FL and FR measurements are fit:
            
        azimark_index = 0            

        if is_FL_FR_refmark_readings_fit( closest_measurements["FL"]["closest-ref-angle"], closest_measurements["FR"]["closest-ref-angle"]  ) :
            # Closest FL/FR measurements fit togehter!
            for FL_FR in ["FL","FR"]:
                key = closest_measurements[FL_FR]["closest-side-key"] 
                closest_index = closest_measurements[FL_FR][key]
                assert not np.isnan( load_data[FL_FR][closest_index]["hctick"] ), "{}".format(load_data[FL_FR][closest_index]["datetime"])
                hc, vc = extract_meas_angle(FL_FR, closest_index)
                result.use_measurement(FL_FR, azimark_index, load_data[FL_FR][closest_index]['datetime'], hc, vc )
        else:
            # Closest FL and FR measurements DO NOT fit togehter! Trying to find a fitted pair, assuming the closest from FL and FR is the refmark reading of this absolute measurement
            del FL_FR
            closest_from_FL_FR = "FL" if closest_measurements["FL"]["diff"] < closest_measurements["FR"]["diff"] else "FR"
            key = closest_measurements[closest_from_FL_FR]["closest-side-key"] 
            closest_index = closest_measurements[closest_from_FL_FR][key]
            hc, vc = extract_meas_angle(closest_from_FL_FR, closest_index)
            result.use_measurement(closest_from_FL_FR, azimark_index, load_data[closest_from_FL_FR][closest_index]['datetime'], hc, vc )
            other_FL_FR = "FL" if closest_from_FL_FR == "FR" else "FL"
            key = "after-index" if closest_measurements[other_FL_FR]["closest-side-key"] == "before-index" else "before-index"
            closest_index = closest_measurements[other_FL_FR][key]
            assert closest_index is not None
            angle_key = "before-ref-angle" if key == "before-index" else "after-ref-angle"
            assert is_FL_FR_refmark_readings_fit( closest_measurements[closest_from_FL_FR]["closest-ref-angle"], closest_measurements[other_FL_FR][angle_key]  ) 
            hc, vc = extract_meas_angle(other_FL_FR, closest_index)
            result.use_measurement(other_FL_FR, azimark_index, load_data[other_FL_FR][closest_index]['datetime'], hc, vc )
            

        result.end_of_data_feed("sep-file")
        # result.print_out()
        return result, pillar_index
        
        

    def check_availability(self):
        assert (not np.isnan(self.last_measured_refmark_true_angle["hc"]["FL"]))  and  (not np.isnan(self.last_measured_refmark_true_angle["hc"]["FR"])), "Reference positions WGS84 data is not supplied"
        assert (not np.isnan(self.last_measured_refmark_theo_angle["hc"]["FL"]))  and (not np.isnan(self.last_measured_refmark_theo_angle["hc"]["FR"])), "Reference position reading is not available"
        
    # Returns the actual (or: last used) average SU/SD angles:
    def act_measured_angles_avg(self):
        
        if self.last_measured_refmark_theo_angle["hc"][self.reference_telescope_position_complement] > self.last_measured_refmark_theo_angle["hc"][self.reference_telescope_position]:
            result = (self.last_measured_refmark_theo_angle["hc"][self.reference_telescope_position] + self.last_measured_refmark_theo_angle["hc"][self.reference_telescope_position_complement] - math.pi) / 2.0
        else:
            result = (self.last_measured_refmark_theo_angle["hc"][self.reference_telescope_position] + self.last_measured_refmark_theo_angle["hc"][self.reference_telescope_position_complement] + math.pi) / 2.0
        
        assert result >= 0.0 and result <=math.tau
        return result
    
    def get_last_measured_vc0(self):
        if (self.last_measured_refmark_true_angle["vc"]["FL"] is not None) and (not np.isnan(self.last_measured_refmark_true_angle["vc"]["FL"])) \
           and \
           (self.last_measured_refmark_theo_angle["vc"]["FL"] is not None) and (not np.isnan(self.last_measured_refmark_theo_angle["vc"]["FL"])):
               vc0_FL = self.last_measured_refmark_true_angle["vc"]["FL"] - self.last_measured_refmark_theo_angle["vc"]["FL"]
               vc0_FR = self.last_measured_refmark_true_angle["vc"]["FR"] - self.last_measured_refmark_theo_angle["vc"]["FR"]
               return np.mean( [vc0_FL, vc0_FR] )
        return None

    # "hpa_t = mod(hpa_n - azimuth_n + azimuth_t, 2*pi)"
    # 'hpa_nt()' does the same
    def h_ned_to_theo(self,h_ned):
        self.check_availability()
        FL_part = (self.last_measured_refmark_theo_angle["hc"]["FL"] - self.last_measured_refmark_true_angle["hc"]["FL"] + math.tau) % math.tau
        FR_part = (self.last_measured_refmark_theo_angle["hc"]["FR"] - self.last_measured_refmark_true_angle["hc"]["FR"] + math.tau) % math.tau
        # SD and SU difference refers to the same WGS-direction/hc-direction combination
        assert abs(FL_part-FR_part)<0.1
        res = h_ned +  angle_utils.angle_avg(FL_part,  FR_part)
        res = (res + math.tau)  % math.tau
        return res

    def h_theo_to_ned(self,h_theo):
        self.check_availability()
        result = h_theo_to_ned(h_theo, self.last_measured_refmark_theo_angle["hc"]["FL"], self.last_measured_refmark_theo_angle["hc"]["FR"], self.last_measured_refmark_true_angle["hc"]["FL"], self.last_measured_refmark_true_angle["hc"]["FR"], False)
        if isinstance(result, tuple): 
            # 'result' is a tuple  of error-message / error-index
            errmsg, error_index = result
            for FL_FR in ["FL","FR"]:
                errmsg += "\n {}: data-source:{} value:{:.2f}deg  date:{}".format(FL_FR, self.meas_data_source[FL_FR], math.degrees(self.last_measured_refmark_theo_angle["hc"][FL_FR]),  mu.exceltime_2_datetimestr(self.meas_data_timestamp[FL_FR]))
            raise ValueError(errmsg)
        return result
    

    # def h_theo_to_ned_from_new_SU(self,h_theo):
    #     if not self.new["SU"]: return math.nan
    #     assert (not np.isnan(self.last_measured_refmark_true_angle["hc"]["SU"]))  and (not np.isnan(self.last_measured_refmark_theo_angle["hc"]["SU"]))
    #     return (h_theo - (self.last_measured_refmark_theo_angle["hc"]["SU"] - self.last_measured_refmark_true_angle["hc"]["SU"]) + mu.pi_x_2) % mu.pi_x_2

    # def h_theo_to_ned_from_new_SD(self,h_theo):
    #     if not self.new["SD"]: return math.nan
    #     assert (not np.isnan(self.last_measured_refmark_true_angle["hc"]["SD"]))  and (not np.isnan(self.last_measured_refmark_theo_angle["hc"]["SD"]))
    #     return (h_theo - (self.last_measured_refmark_theo_angle["hc"]["SD"] - self.last_measured_refmark_true_angle["hc"]["SD"]) + mu.pi_x_2) % mu.pi_x_2

    # Returns the approximate position of the NED North direction.
    # in theo reference frame, when vc is the given angle.
    # def hc_theo_SU_north(self,theo_vc, instr):
    #     assert self.reference_sensor_position == instr.hc_angle_reference_sensor_pos()
    #     # Note: 'instr' parameter can be deleted!

    #     assert not np.isnan(self.last_measured_refmark_theo_angle["hc"]["SU"])
    #     assert not np.isnan(self.last_measured_refmark_true_angle["hc"]["SU"])
    #     result = mu.pi_x_2 - (self.last_measured_refmark_true_angle["hc"]["SU"] - self.last_measured_refmark_theo_angle["hc"]["SU"])
    #     if not instr.is_sensor_up_at_vc(theo_vc):
    #     # if (theo_vc + mu.pi_x_2) % mu.pi_x_2 > np.pi:  # Theo is in 180...360 deg range
    #         result = result + np.pi
    #     result = (result + mu.pi_x_2) % mu.pi_x_2
    #     return result



# Converts THEO vertical angle to NED vertical angle.
# 'vc_NED=0° means horizontal.
# 'vc_THEO=0°' means vertical, objective pointing upward.
def v_ned_to_theo(v_ned):
    return (v_ned + math.pi/2) % math.tau
    # return mu.pi_x_2 - (v_ned + math.pi/2) % mu.pi_x_2

def v_theo_to_ned(v_theo):
    return (v_theo - math.pi/2 + math.tau) % math.tau


# 'default_azimuth_def' must be set for LM callback function, that can not receive a class instance through function arguments
default_azimuth_def = None
def set_default_azimuth_def(azi_def):
    global default_azimuth_def
    default_azimuth_def = azi_def


# Reads in reference mark readings from DAT files. Used to read ADS reference mark data files.
def read_refmark_readings_from_dat_files(data_source_def, pillar_index, azimark_index, for_timelabels, meas_pos, instr, COMPNAME__PILLAR_INDEX):
    assert isinstance(data_source_def, azimark_measured_data_source_def_class)
    assert data_source_def.source_type == azimark_meas_data_source_type_enum.separate_dat_files
    
    assert len(meas_pos.pillars) == 1, "This function features only 1 pillar"
    
    first_full_filename = None
    load_dat_file_name_list = {  "FL": [], "FR": [] }
    load_dat = {  "FL": None, "FR": None }

    def load_dat_file(for_timelabel, must_exist):
        nonlocal load_dat_file_name_list        
        found_files = 0
        for SU_SD in ["SU","SD"]:
            FL_FR = instr.SUSD_2_FLFR(SU_SD)
            filename = data_source_def.file_name_def
            datetime = mu.exceltime_2_datetime(for_timelabel)
            filename = datetime.strftime(filename)
            # ToDo: make the replacement only once in the parent function
            filename = gu.replace_masks_in_filename(filename, meas_pos, pillar_index) # Replacing pier name
            filename = filename.replace(gu.MASK_SENSORPOS, SU_SD)
            filename = filename.replace(gu.MASK_THEO_FACE_POSITION, FL_FR )
            assert data_source_def.azimark_name is not None
            filename = filename.replace(gu.MASK_REFPOS, data_source_def.azimark_name)
            
            nonlocal first_full_filename
            if first_full_filename is None:
                first_full_filename = filename
                # Checking if reference position data file directory exist:
                data_dir = os.path.dirname(filename)
                assert os.path.isdir(data_dir), "Azimuth data file directory \"{}\" does not exist".format(data_dir)
                

            if filename in load_dat_file_name_list[FL_FR]: 
                found_files += 1
                continue

            filename_only = os.path.basename(filename)
            
            # ToDo: optimize file searching            
            if not os.path.isfile(filename):
                assert not must_exist, "Reference mark reading dat file \"{:}\" not found".format(filename)
                # print(" no such file")
                continue

            found_files += 1
            load_dat_file_name_list[FL_FR] += [filename]
            print(f"   {filename_only} : ", end = '')

            d = np.genfromtxt(filename, names=True , usecols=[0,5,6], ndmin = 1)
            
            d = d[  ~np.isnan(d['hctick'])  ]  # Removing records without results
            if d.shape[0]==0:
                print("file has no valid record")
                continue

            print("{} {} read".format( d.shape[0], mu.plural_suffix(d.shape[0],"record","records")    ))
            assert "datetime" in d.dtype.names
            assert "hctick" in d.dtype.names
            assert "vctick" in d.dtype.names
            if load_dat[FL_FR] is None:
                load_dat[FL_FR] = d # First of this type of data
            else:
                load_dat[FL_FR] = np.concatenate( (load_dat[FL_FR],d) )
                load_dat[FL_FR].sort(order="datetime")                
        return found_files > 0
                
    
    timelabel_from = for_timelabels.min()
    timelabel_to   = for_timelabels.max()
    
    mask_datetime_interval = mu.mask_to_date_interval(data_source_def.file_name_def)
    date_from = timelabel_from - data_source_def.max_time_diff_to_use_a_refmark_reading__days
    date_to   = timelabel_to   + data_source_def.max_time_diff_to_use_a_refmark_reading__days
    filename_dates_to_load = mu.generate_date_range(date_from, date_to, mask_datetime_interval )
    assert len(filename_dates_to_load) >= 1
    
    load_files = 0
    date_range_str = "{}...{}".format( mu.exceltime_2_datetimestr(date_from, geomag_utils_preferences.UI_DATETIME_FORMAT_MIN), mu.exceltime_2_datetimestr(date_to, geomag_utils_preferences.UI_DATETIME_FORMAT_MIN)  )    
    print("Loading azimuth mark data files from directory \"{}\" for date range of {}:".format(os.path.dirname(data_source_def.file_name_def),date_range_str))

    for filename_date in filename_dates_to_load:
        if load_dat_file(filename_date, False):
            load_files += 1
    assert (load_files  > 0) , "No reference data found for the date range \"{}\" (searched file name example:\"{}\")".format(date_range_str,first_full_filename) 

    missing_data_positions = ""
    if (load_dat['FL'] is None): missing_data_positions = "FL"
    if (load_dat['FR'] is None): missing_data_positions = "FR" if missing_data_positions == "" else "FL+FR"
    assert  missing_data_positions == "", "No {} reference mark data found for date range of {}".format(missing_data_positions,date_range_str)

    print("{} FL and {} FL azimuth data files loaded".format(  len(load_dat_file_name_list["FL"]), len(load_dat_file_name_list["FR"])) )
    
    for FL_FR in ["FL", "FR"]:
         load_dat[FL_FR] = rcf.append_fields(load_dat[FL_FR], (COMPNAME__PILLAR_INDEX, COMPNAME__AZIMARK_INDEX), ([0],[0]),  (np.uint32,np.uint32), usemask=False )
         load_dat[FL_FR][COMPNAME__PILLAR_INDEX][:] = pillar_index
         load_dat[FL_FR][COMPNAME__AZIMARK_INDEX][:] = azimark_index
         
    
    return load_dat
    
# default_azimuth_def               = azimuth_def_class("internal-table" , [ "dat-file","external-table"] )  #
# azimuth_def_without_dat_file_data = azimuth_def_class("internal-table" , [ "external-table"] )  #




