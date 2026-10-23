# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 08:24:45 2024

@author: Merényi László

Imports absolute geomagnetic measurements data from XLS(X) files
Currently supported files:
- Mingeo manual filled XLSX
- Mingeo DS XLSX
- CYG XLSX

"""


import pandas as pd
import datetime
import numpy as np
import math
import sys
import glob
import os
import tabulate
import fnmatch
from enum import Enum, auto

import misc_utils as mu
import angle_utils
import amdp_base_formulas as amdp
import geomag_utils as gu
import dim_utils
import dim_measurement_loader 
import cli_utils as cu


NR_OF_DIM_MEAS_MAIN_POSITIONS = 9    # 4x declination + 4x inclination + scale-check
NR_OF_DIM_REFMARK_MEAS_POSITIONS = 4  # FL/FR at start and end
REFMARK_READING_ESTIMATED_DURATION__SEC = 20  # Used for time stamp calculation if reference mark reading time is not noted

def_measurement_date_is_read_from_cell = False    # If 'True' then XLS files that seems to be out of date limit according to their file name are still read-in,
                                              # and date is read from XLS cell, then finel date check is done by the cell date.
                                              # This results is prolonged processing time, but this mode is useful to check filename and cell date's consistency.
                                              #

supported_input_file_types = ("Mingeo-manual-XLSX", "Mingeo-DS-XLSX", "CYG-XLSX")

class DI_spreadsheet_file_type_enum(Enum):
    Mingeo_manual_XLSX = supported_input_file_types[0]
    Mingeo_DS_XLSX     = supported_input_file_types[1]
    CYG_XLSX           = supported_input_file_types[2]
    

PILLAR_IN_FILENAME_MASK  = "[%PILLAR]"
INSTR_IN_FILENAME_MASK   = "[%INSTR]"
STATION_IN_FILENAME_MASK = "[%STATION]"

# Some main definition of the different file formats, describing file name format and worksheet (tabs) usaged
input_file_defs =  \
{  
   "Mingeo-manual-XLSX": {  "default-filename-mask" : STATION_IN_FILENAME_MASK + "*_%Y-%m-%d_"+INSTR_IN_FILENAME_MASK+".xlsx", "pillar-subdir":"", "worksheet_mask" : "my sequence*" },  # AFter "MNY_" theree can be optionally the pillar
   "Mingeo-DS-XLSX":     {  "default-filename-mask" : INSTR_IN_FILENAME_MASK + "_%Y%m%d.xlsx", "pillar-subdir":"", "worksheet_mask" : "*"  } ,  # In some files: "*_%Y%m%d"
   "CYG-XLSX":           {  "default-filename_mask" : "CYG baseline caculation_%Y.%m.%d*.xlsx", "pillar-subdir":"", "worksheet_mask" : "CYG%Y%m%d_*" }
}

# Detailed descriptions of cells data
# XLSX cell's reference: 
#  - first index:  row number, starting from 0, so excel 'Row 1' is 0 
#  - second index: column number, 0 meaning Excel/Calc "A" column
xlsx_cells_defs = \
{
    DI_spreadsheet_file_type_enum.Mingeo_manual_XLSX: {   # XLSX format of manual measuremnets
        
        "checks": [ [37,0,'Inclination series readings'] ],
        
        "date"     : [1,2],
        "timezone" : [2,2],
        "IMO-code" : [4,2],
        "refmark_descr" : [6,2],
        "refmark-WGS84" : [9,10],
        "refmark-WGS84_unit" : "deg",
        "pillar"   : [12,2],
        "theo_SN" : [18,2],
        "electronics_SN" : [22,2],
        "angleunit": [24,2],
        "decl"     : [36,16],
        "incl"     : [43,17],
        "notused_starttime": [32,0], 
        "notused_endtime": [42,0] ,
        
        "NM_D_measpart_row_index" : 32,
        "NM_I_measpart_row_index" : 38,
        "NM_SC_measpart_row_index" : 42,
        
        "NM_measpart_rotated_axis_col_index" : None, 
        "NM_measpart_timestamp_col_index" : 0 ,
        "NM_measpart_M_col_index"         : 10,
        
        "start-MIRA-SU"                   : 29,
        "start-MIRA-SD"                   : 30,
        "end-MIRA-SU"                     : 45,
        "end-MIRA-SD"                     : 46,
        
        "NM_measpart_posname_col_index"   : None

        },

    DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX: {  # DS instrument saved XLSX format 
        
        "timezone"       : "UTC",        
        "angleunit":       "dms",
        "checks": [ [0,0,'Measurement ID:'] ],
        
        "date"           : [1,1, "Date:"],
        "theo_SN"        : [7,1, 'Theodolite serial:'],
        "electronics_SN" : [11,1, 'Electronics serial:'],

        "IMO-code"         : [21,1, 'IMO name/code:'],
        "pillar"           : [23,1, 'Pillar name/code:'],        
        "refmark_descr"    : [25,1, 'Azimuth mark name:'],
        "refmark-WGS84"    : [26,1, 'Azimuth mark angle:'],

        "decl"             : [48,1, "Declination, D:"],
        "incl"             : [49,1, "Inclination, I:"],
        
        "NM_D_measpart_row_index" : 39,
        "NM_I_measpart_row_index" : 35,
        "NM_SC_measpart_row_index" : 34,


        "NM_measpart_rotated_axis_col_index" : None, 
        "NM_measpart_timestamp_col_index" : 6,
        "NM_measpart_M_col_index"         : 4,
        
        "start-MIRA-SD"                   : [43,2,"MD"],
        "start-MIRA-SU"                   : [44,2,"MU"],
        "end-MIRA-SD"                     : [32,2,"MD"],
        "end-MIRA-SU"                     : [33,2,"MU"],
        
        "NM_measpart_posname_col_index"   : 1
        
        },
    
    DI_spreadsheet_file_type_enum.CYG_XLSX:  { # XLSX files used by Saevit
        "timezone"    : "UTC",        
        "angleunit"   : "gon",
        
        "checks": [ [1,0,'NULL METHOD'], [3,1,'azimuth'] ],
        
        "start-MIRA-SD"                   : [8,1,"Mark Down"],
        "start-MIRA-SU"                   : [7,1,"Mark Up"],
        "end-MIRA-SD"                     : [10,1,"Mark Down"],
        "end-MIRA-SU"                     : [9,1,"Mark Up"],
        
        "IMO-code"         :     "CYG",
        "electronics_SN"   :     "xxx",
        "theo_SN"          :     "xxx",
        "pillar"           :     "%obs-default", 
        "refmark_descr"    :     "xxx",
        "refmark-WGS84"    :     [ 4,1],
        "decl"             : [19,1],
        "incl"             : [28,1],
        
        "NM_D_measpart_row_index"  : 15,
        "NM_I_measpart_row_index"  : 24,
        "NM_SC_measpart_row_index" : None,

        "NM_measpart_rotated_axis_col_index" : 1, 
        "NM_measpart_timestamp_col_index" : 2,
        "NM_measpart_M_col_index"         : None,
        
        "NM_measpart_posname_col_index"   : 0
        
    }
 
}
Mingeo_XLSX_pos_order = ("ED", "WD", "EU", "WU", "NU", "SD", "ND", "SU", "SC")


def check_angle_unit_str(angle_unit_str):
    
    if angle_unit_str == "degmin" : angle_unit_str = "dm"
    if angle_unit_str == "d": angle_unit_str = "deg"
    assert  angle_unit_str in  ["gon", "deg", "dms", "dm" ], "Unsupported angle unit: {}".format(angle_unit_str)
    
    return angle_unit_str



def read_absmeas_XLS_file(input_file_type, XLS_file_name, worksheet_name, meas_pos, date_from_filename, theo_SN_filter = None, datetime_filter = None, DIM_instr = None, consistency_checks = True): 
    """ Reads in an and processes one absolute measurement from an XLS(X) file.
    
        If an XLS(X) file has more worksheets with separate measurements, then it reads only one from them
        'input_file_type' tells the format of the XLSX file.
    """
    assert input_file_type in supported_input_file_types
    input_file_type = DI_spreadsheet_file_type_enum(input_file_type)
    xlsx_cells_def = xlsx_cells_defs[input_file_type]
    
    def check_cell( check_res, cell_coordinates_def, error_msg):
        if check_res: return
        if isinstance (cell_coordinates_def, str) : 
            cell_coordinates = xlsx_cells_def[cell_coordinates_def]
        else:
            assert len(cell_coordinates_def) in [2,3], cell_coordinates_def  # Normally 'cell_coordinates_def' has a pair of coordinates, but "check" cells definition has a thrird one: the expected value
            cell_coordinates = cell_coordinates_def[:2]
            
        error_msg = error_msg+"! File:{} Cell:{}:{}".format( os.path.basename(XLS_file_name), chr(ord('A')+cell_coordinates[1]), cell_coordinates[0]+1) 
        os.startfile(XLS_file_name)
        assert False, error_msg
    
    def read_cell(cell_id, force_result_type = None, can_be_empty = False):
        nonlocal xlsx_cells_def
        cell_def = xlsx_cells_def[cell_id]
        if isinstance(cell_def,list):
            # 'cell_def' refers to an XLSX cell
            row_nr = cell_def[0]
            col_nr = cell_def[1]
            result = absmeas_table.loc[ row_nr, col_nr]
            if len(cell_def)>2:
                check_cell_label = cell_def[2]
                assert check_cell_label == absmeas_table.loc[ row_nr, col_nr-1]
        else:
            # 'cell_def' is a constant value - so not read from XLSX
            result = cell_def
        if isinstance(result,float) and np.isnan(result):
            if can_be_empty: return None
            check_cell (False, cell_def, "Cell 'is empty")
        if (force_result_type is not None):
            if force_result_type == "str":
                if not isinstance(result, str):
                    result = "{}".format(result)
            elif force_result_type == "float":
                if isinstance(result, str):
                    result = float(result)
                else:
                    assert isinstance(result, float)
            else:
                assert False
            
        return result
    
    # Converts a timestamp cell containing 'hour-minute' or 'hour-minute-second' data into full datetime.
    # It also converts result into UTC time
    def cell_timestamp_to_datetime(cell_pos, use_datetime_if_empty = None ):
        nonlocal timezone_corr_hours, measurement_day
        timestamp = absmeas_table.loc[ cell_pos[0] ][  cell_pos[1]  ]
        if isinstance(timestamp,float) and np.isnan(timestamp):
            check_cell( (use_datetime_if_empty is not None) and (isinstance(use_datetime_if_empty, datetime.datetime)), cell_pos, "DateTime is not available")
            
            return use_datetime_if_empty
        
        if isinstance(timestamp, str) and (":" in timestamp):
            # "hh:mm:ss" format is used in DS files
            hours = int(timestamp[:timestamp.index(":")])
            _timestamp = timestamp[timestamp.index(":")+1:]
            minutes = int(_timestamp[:_timestamp.index(":")])
            _timestamp = _timestamp[_timestamp.index(":")+1:]
            seconds = int(_timestamp)
        else:
            # Manual measurement timestamp is formatted as "hhmmss"
            if isinstance(timestamp, str): timestamp = int(timestamp)            
            assert timestamp>=0
            seconds = 0
            if timestamp > 2359:
                assert timestamp > 9999
                seconds = timestamp % 100
                timestamp = timestamp // 100
            minutes = timestamp % 100
            timestamp = timestamp // 100
            hours = timestamp
        check_cell( (seconds<60) and (minutes<60) and (hours<24), cell_pos, "Wrong time: \"{}\"".format(timestamp))

        measurement_time = datetime.time(hours, minutes, seconds)

        result_datetime = datetime.datetime.combine(measurement_day, measurement_time)
        result_datetime = result_datetime -  datetime.timedelta(hours=timezone_corr_hours)
        
        return result_datetime
    
    def convert_to_deg(angle) :
        nonlocal angle_unit_str
        if angle_unit_str == "deg": 
            return angle
        elif angle_unit_str == "gon":
            return angle/400.0*360.0
        elif angle_unit_str == "dms":
            assert type(angle) == str
            deg_pos = angle.index("°")
            min_pos = angle.index("'")
            sec_pos = angle.index("\"")
            assert sec_pos > min_pos + 2 > deg_pos + 2
            degree   = int(angle[:deg_pos])
            minute   = int(angle[deg_pos+1:min_pos])
            second   = float(angle[min_pos+1:sec_pos]) # Second can have decimal digit
            assert 0 < minute < 60
            assert 0 < second < 60
            return degree + minute / 60.0 + second / 3600.0
            
        else:
            assert False, "not yet implemented {} unit conversion for \"{}\"".format(angle_unit_str, angle)
            
    def convert_to_rad(angle) :
        nonlocal angle_unit_str
        if angle_unit_str == "deg": 
            return math.radians(angle)
        elif  angle_unit_str == "gon":
            return math.radians(angle/400.0*360.0)
        else:
            assert False, "not yet implemented for {}".format(angle_unit_str)
                
    def DS_xlsx_file__convert_to_deg( cell_str ):
        # 2026.01.26: for general data, the format is always "<d>° <m>' <s>"", while for raw data, format can be the
        #   "<d>° <m>' <s>"", or "gon" or "<n>°"
        nonlocal angle_unit_str
        
        if not '°' in cell_str:
            # "gg.gggg" gon format
            assert angle_unit_str == "gon", f"${angle_unit_str=} ${cell_str=}"
            assert "." in cell_str
            assert not "\'" in cell_str
            return float(cell_str) / 400.0 * 360.0 

        if cell_str[-1] == "°":
            degrees = float(cell_str[:-1])
            return degrees
        
        else:
            assert angle_unit_str == "dms"
            
            if not cell_str[-1] in [ '"', '”' ]: # Last char should be the second unit
                assert False, cell_str
            deg_unit_pos = cell_str.index("°")
            min_unit_pos = cell_str.index("\'") 
            degree_str =  cell_str[ :deg_unit_pos ]
            minute_str =  cell_str[ deg_unit_pos+2:min_unit_pos ]
            minutes  = int(minute_str)
            assert minutes>=0 and minutes<60
            second_str =  cell_str[ min_unit_pos+2 : -1 ]
            seconds = float(second_str)
            assert seconds >=0 and seconds < 60.0
            
            degrees = int(degree_str)
            sign = -1 if degrees < 0.0 else +1
            return degrees + sign * minutes/60.0 + sign * seconds / 3600.0
        

    def manual_xlsx_file__parse_angles(row_index, series_part):
        nonlocal angle_unit_str
        if series_part in [ "Decl", "MIRA" ]:
            rotated_angle_degree_input_col_index = 2
            fixed_angle_degree_input_col_index  = rotated_angle_degree_input_col_index + 4
            rotated_angle_col_index__always_deg = 16
            fixed_angle_col_index__always_deg  = 17
        elif series_part == "Incl":
            rotated_angle_degree_input_col_index = 6
            fixed_angle_degree_input_col_index  = rotated_angle_degree_input_col_index - 4
            rotated_angle_col_index__always_deg = 17
            fixed_angle_col_index__always_deg = 16
        else:
            assert False

        if angle_unit_str == "gon":
            rotated_angle_deg = float(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index+3]) / 400.0 * 360.0
            fixed_angle_deg   = float(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+3]) / 400.0 * 360.0
        elif angle_unit_str == "deg":
            rotated_angle_deg = float(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index]  )
            fixed_angle_deg = float(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index]  )
        elif angle_unit_str == "dm":
            _deg = float(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index+0])
            _min = float(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index+1])
            assert 0 <= _deg < 360
            assert 0 <= _min < 60
            rotated_angle_deg = _deg + _min/60.0
            _deg = float(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+0])
            _min = float(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+1])
            assert 0 <= _deg < 360
            assert 0 <= _min < 60, "min={}".format(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+1])
            fixed_angle_deg = _deg + _min/60.0
        elif angle_unit_str == "dms":
            _deg = int(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index+0])
            _min = int(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index+1])
            _sec = int(absmeas_table.loc[row_index, rotated_angle_degree_input_col_index+2])
            assert 0 <= _deg < 360
            assert 0 <= _min < 60
            assert 0 <= _sec < 60
            rotated_angle_deg = _deg + _min/60.0 + _sec/3600.0
            _deg = int(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+0])
            _min = int(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+1])
            _sec = int(absmeas_table.loc[row_index, fixed_angle_degree_input_col_index+2])
            assert 0 <= _deg < 360
            assert 0 <= _min < 60
            assert 0 <= _sec < 60
            fixed_angle_deg = _deg + _min/60.0 + _sec/3600.0
        else:
            assert False

        # Other, more easy way is to use the "always degree" column
        rotated_angle_excel_calculated_deg = float(absmeas_table.loc[row_index, rotated_angle_col_index__always_deg])
        # ---------------------------------------------------------------------------------------------------------------------------------------
        # Note: XLSX with residual correction (from May 2025), the check angle can contain the residual correction, and it can be different than read angle!!!!
        # Take care of this situation!
        # ---------------------------------------------------------------------------------------------------------------------------------------

        # if not abs(rotated_angle_excel_calculated_deg - rotated_angle_deg)<0.01:
        #    assert False, "{:}: {:.3f}  {:.3f}".format(XLS_file_name,rotated_angle_excel_calculated_deg,rotated_angle_deg)
        # rotated_angle_deg = rotated_angle_excel_calculated_deg
        
        fixed_angle_excel_calculated_deg = float(absmeas_table.loc[row_index, fixed_angle_col_index__always_deg])
        assert abs(fixed_angle_excel_calculated_deg - fixed_angle_deg)<1.0
        fixed_angle_deg = fixed_angle_excel_calculated_deg

        assert not np.isnan(rotated_angle_deg), "{}".format(rotated_angle_deg)
        return math.radians(rotated_angle_deg), math.radians(fixed_angle_deg)
        
    
    def DS_xlsx_file__parse_NMpart_angles(row_index, series_part):
        hc_angle = absmeas_table.loc[row_index, 2]
        vc_angle = absmeas_table.loc[row_index, 3]
        hc_angle_deg = DS_xlsx_file__convert_to_deg(hc_angle)
        vc_angle_deg = DS_xlsx_file__convert_to_deg(vc_angle)
        if series_part in [ "Decl", "MIRA" ]:
            rotated_angle_deg = hc_angle_deg
            fixed_angle_deg   = vc_angle_deg
        elif series_part == "Incl":
            rotated_angle_deg = vc_angle_deg
            fixed_angle_deg   = hc_angle_deg
        else:
            assert False
        return math.radians(rotated_angle_deg), math.radians(fixed_angle_deg)
    
    # Read rotated angle from a cell, in a more generaly form
    # If fixed axis is not written in Excel, then this function determines it.
    def gen_parse_NMpart_angles(series_part, pos_name, rotated_axis_row_index, rotated_axis_col_index, fixed_axis_col_index = None, D_from_decl_series = None):
        # Developed for "CYG" format in 2026. ToDo: extend it ...
        nonlocal azimuth_mark_NED, MIRA_mean_for_SU_pos
        rotated_angle = absmeas_table.loc[rotated_axis_row_index, rotated_axis_col_index]
        rotated_angle = convert_to_rad(rotated_angle)
        if fixed_axis_col_index is None:
            # Fixed axis is not written into the talbe
            if series_part == "Decl":
                if pos_name[1]=='U':
                    fixed_angle = np.pi / 2.
                elif pos_name[1]=='D':
                    fixed_angle = 3 * np.pi / 2.
                else: 
                    assert False
            elif series_part == "Incl":
                assert D_from_decl_series is not None
                D_theo = amdp.hpa_nt(D_from_decl_series, azimuth_mark_NED, MIRA_mean_for_SU_pos  ) 
                if pos_name in [ "NU","SD"]:
                    fixed_angle =  D_theo
                elif pos_name in [ "ND", "SU"]:
                    fixed_angle = (D_theo + np.pi) % (2*np.pi)
                else:
                    assert False
                
            else:
                assert False
        else:
            assert False,"Not implemented"
            
        return rotated_angle, fixed_angle
        
        
    absmeas_table = pd.read_excel(XLS_file_name, worksheet_name, header=None, index_col=None )
    
    # Starting values of basic solution:
    delta_declination_part = 0.0
    epsilon = 0.0
    M_offset_declination_part = 0.0
    M_scale = 1.0
    
    for checking_cell in xlsx_cells_def['checks']:
        check_cell( absmeas_table.loc[ checking_cell[0] ] [ checking_cell[1] ] == checking_cell[2], checking_cell, "Excel cell check error: CELL[{}:{}]={} <> {}".format( checking_cell[0], checking_cell[1], absmeas_table.loc[ checking_cell[0] ] [ checking_cell[1] ] , checking_cell[2]) )
    
    IMO_code_from_XLS = read_cell( "IMO-code" ) 
    IMO_code_from_XLS = IMO_code_from_XLS.upper()
    
    assert  IMO_code_from_XLS == meas_pos.station_code, "IMO-code-from-XLS:{}  station:{}".format(IMO_code_from_XLS, meas_pos.station_code)  # E.g. "MNY"

    angle_unit_str = read_cell('angleunit')  # In Mingeo-manual-XLSX, the same unit is used in all data
                                             # In DS, different units may be used in case of general data and raw measurement and azimuth mark dir.
                                             
    angle_unit_str = check_angle_unit_str(angle_unit_str)


    electronics_SN = read_cell('electronics_SN', "str") 
    theo_SN        = read_cell('theo_SN', "str" )

    if (theo_SN_filter is not None) and (not theo_SN in theo_SN_filter):
        return "Theo SN \"{}\" is not in accepted SN list".format(theo_SN), [None,None,None,None],None,[None,None],None,None,[None,None,None,None],None

    timezone_str = read_cell("timezone")
    timezone_str = timezone_str.upper()
    timezones = { "UTC":0, "CET":1, "CEST":2, "KST": 9 }
    
    assert timezone_str in timezones.keys(), "Unsupported time-zone: {}".format(timezone_str)
    timezone_corr_hours = timezones[timezone_str]
   
    
    if xlsx_cells_def['pillar'] == "%obs-default":
        pillar = gu.measurement_site_list[meas_pos.station_code]["piers"][0]
    else:
        pillar =   read_cell('pillar', None)
    
    refmark_descr = read_cell('refmark_descr')
    
    refmark_WGS84 = read_cell('refmark-WGS84')
    if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX:
        azimuth_mark_NED_deg = DS_xlsx_file__convert_to_deg(refmark_WGS84)
    elif "refmark-WGS84_unit" in xlsx_cells_def.keys():
        assert xlsx_cells_def["refmark-WGS84_unit"] == "deg" # Mingeo XLSX file contains refmark WGS84 in selected angle and in deg
        azimuth_mark_NED_deg = refmark_WGS84
    else:
        azimuth_mark_NED_deg = convert_to_deg(refmark_WGS84) 
    if isinstance(azimuth_mark_NED_deg, str):
        azimuth_mark_NED_deg=float(azimuth_mark_NED_deg)
    azimuth_mark_NED = math.radians(azimuth_mark_NED_deg)
    
    assert isinstance(azimuth_mark_NED_deg,float), azimuth_mark_NED_deg
    if meas_pos.station_code == "MNY":
        if refmark_descr[:11] == "Water Tower": refmark_descr = "TowerLeft" # Short name of water tower

    if "date" in xlsx_cells_def.keys():
        measurement_day = read_cell("date")
    else:
        measurement_day  = date_from_filename  # Date is from filename

    if isinstance(measurement_day, str):
        # Sometimes the date cell is read out as string, sometime as datetime.datetime class
        measurement_day  = measurement_day.replace("."," ")
        measurement_day  = datetime.datetime.strptime(measurement_day,'%Y %m %d')

    assert isinstance(measurement_day, datetime.datetime), f"${XLS_file_name=} ${measurement_day=}"
    assert (measurement_day.year>2020) and (measurement_day.hour == 0) and (measurement_day.minute == 0)
    # Using only the date information:

    # Checking if file name date corresponds to date stored in the file
    if date_from_filename is not None:
        date_diff = measurement_day - date_from_filename
        assert abs(date_diff.days) <= 1, "File date error in file {:}: {:} - {:}".format(XLS_file_name, measurement_day, date_from_filename )

    measurement_day = measurement_day.date()
    
    # In DS files, I/D may be reported in different angle unit than the raw measurement, so reading out I/D now, before the raw measurements
    # In DS, "I_raw"/"D_raw" cells can be text!
    I_raw = read_cell('incl', "float" if input_file_type != DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX else None)
    
    # In 2026.01 january MNY DS files, the D values are wrong, and are deleted from cells, so we allow empty (None) values for DS files
    D_raw = read_cell('decl', "float" if input_file_type != DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX else None, input_file_type == DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX)
    
    if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX:
        # In DS, I and D are given in three angle units. Ensuring, that we are using the good one
        row_nr = xlsx_cells_def['decl'][0]
        col_nr = xlsx_cells_def['decl'][1]
        assert absmeas_table.loc[ row_nr-1, col_nr] == "dms"
        I_from_file_deg = DS_xlsx_file__convert_to_deg(I_raw)
        D_from_file_deg = None if D_raw is None else DS_xlsx_file__convert_to_deg(D_raw)
    elif input_file_type == DI_spreadsheet_file_type_enum.Mingeo_manual_XLSX:
        I_from_file_deg, D_from_file_deg = I_raw, D_raw
    else:
        I_from_file_deg, D_from_file_deg = convert_to_deg(I_raw), convert_to_deg(D_raw)
    I_from_file = math.radians(I_from_file_deg)        
    D_from_file = None if D_from_file_deg is None else math.radians(D_from_file_deg) 
   

    # Reading in declination series points (time and angles)
    NM_D_measpart_first_row_index = xlsx_cells_def["NM_D_measpart_row_index"] 
    NM_I_measpart_first_row_index = xlsx_cells_def["NM_I_measpart_row_index"] 
    NM_SC_measpart_row_index = xlsx_cells_def["NM_SC_measpart_row_index"] 
    
    NM_measpart_act_row_index = NM_D_measpart_first_row_index
    
    NM_measpart_rotated_axis_col_index = xlsx_cells_def["NM_measpart_rotated_axis_col_index"] 
    NM_measpart_posname_colindex  = xlsx_cells_def["NM_measpart_posname_col_index"] 
    timestamp_input_col_index = xlsx_cells_def["NM_measpart_timestamp_col_index"] 
    M_colindex = xlsx_cells_def["NM_measpart_M_col_index"]  # M is non zero in residual measurements!


    position_dtype =  [  (dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION, dim_measurement_loader.COMPTYPE__TELESCOPE_CARDINAL_DIRECTION), (dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION, dim_measurement_loader.COMPTYPE__TELESCOPE_FACE_POSITION) ]

    mag_meas_positions = np.empty( 9, dtype = position_dtype +  [ ("pos_name", "S2"), ("hc",float),  ("vc",float), ("timelabel_datetime", datetime.datetime), ("timelabel", float), ("M",float) ]   )
    azimark_positions  = np.empty( 9, dtype = position_dtype +  [ ("pos_name", "U2"), ("hc",float),  ("vc",float), ("timelabel_datetime", datetime.datetime), ("timelabel", float) ] )
    
    act_pos_index = 0
    D_series_timelabels = []
    D_series_indeces = [] 
    I_series_timelabels = []
    I_series_indeces = []
    posindex = {}
    
    if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX:
        unit_def = absmeas_table.loc[ 31, 2]  # Should be "H angle (gon)" or "H angle (d m s)"
        assert unit_def[:9]=="H angle ("
        if unit_def[9:] == "gon)":
            angle_unit_str = "gon"
        elif unit_def[9:] == "d m s)":
            angle_unit_str = "dms"
        else:
            assert False
    
    if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_manual_XLSX:    
        start_MIRA_UP, sink = manual_xlsx_file__parse_angles( xlsx_cells_def["start-MIRA-SU"], "MIRA" )
        start_MIRA_DN, sink = manual_xlsx_file__parse_angles( xlsx_cells_def["start-MIRA-SD"], "MIRA" )
        end_MIRA_UP, sink   = manual_xlsx_file__parse_angles( xlsx_cells_def["end-MIRA-SU"], "MIRA" )
        end_MIRA_DN, sink   = manual_xlsx_file__parse_angles( xlsx_cells_def["end-MIRA-SD"], "MIRA" )
    else:    
        start_MIRA_DN    = read_cell( "start-MIRA-SD" )
        start_MIRA_UP    = read_cell( "start-MIRA-SU" )
        end_MIRA_DN      = read_cell( "end-MIRA-SD" )
        end_MIRA_UP      = read_cell( "end-MIRA-SU" )
        if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX:
            end_MIRA_DN   = math.radians(DS_xlsx_file__convert_to_deg( end_MIRA_DN ))        
            start_MIRA_DN = math.radians(DS_xlsx_file__convert_to_deg( start_MIRA_DN ))        
            start_MIRA_UP = math.radians(DS_xlsx_file__convert_to_deg( start_MIRA_UP ))
            end_MIRA_UP   = math.radians(DS_xlsx_file__convert_to_deg( end_MIRA_UP ))        
        else:
            start_MIRA_DN,start_MIRA_UP,end_MIRA_DN,end_MIRA_UP =  convert_to_rad(start_MIRA_DN), convert_to_rad(start_MIRA_UP), convert_to_rad(end_MIRA_DN), convert_to_rad(end_MIRA_UP)
           

    # Nem megy: check_cell( abs(start_MIRA_UP - end_MIRA_UP) < math.radians(0.1), "start-MIRA-SU" , "Start and end MIRA UP angle mismatch" )
    assert abs(start_MIRA_UP - end_MIRA_UP) < math.radians(0.1)
    assert abs(start_MIRA_DN - end_MIRA_DN) < math.radians(0.1), "Start-MIRA-DN:{:.2f}°   End-MIRA-DN:{:.2f}° ".format( math.degrees(start_MIRA_DN), math.degrees(end_MIRA_DN) )
    assert abs(angle_utils.angle_distance(start_MIRA_UP, start_MIRA_DN)) - np.pi < 0.05
    

    def read_in_next_meas_part(series_part, estiamte_timelabel_if_missing = False):
        nonlocal mag_meas_positions, act_pos_index, NM_measpart_act_row_index, NM_measpart_rotated_axis_col_index, D_series_timelabels,D_series_indeces,I_series_timelabels, I_series_indeces,  M_colindex
        assert series_part in ["Decl", "Incl"]
        
        if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_manual_XLSX:
            pos_name = Mingeo_XLSX_pos_order[act_pos_index]
        else:
            # Looking for line of 'act_pos_index'
            # NM_measpart_act_row_index = -1
            #for i in range(9):
            #    xlsx_pos_name =  absmeas_table.loc[NM_measpart_first_row_index + i, 1] 
            #    if xlsx_pos_name == "CH": xlsx_pos_name = "SC"
            #    if xlsx_pos_name == pos_name:
            #        NM_measpart_act_row_index = NM_measpart_first_row_index + i
            #        break
            #assert NM_measpart_act_row_index > 0
            pos_name = absmeas_table.loc[NM_measpart_act_row_index, NM_measpart_posname_colindex] 
            if pos_name == "CH": pos_name = "SC" # In case of DS, scale-check is called "CH"
            if len(pos_name)>2: 
                # "East Up"
                pos_name_parts = pos_name.split(' ')
                assert len(pos_name_parts) == 2
                pos_name = pos_name_parts[0][0] + pos_name_parts[1][0]
            if series_part == "Decl":
                assert pos_name in Mingeo_XLSX_pos_order[0:4]
            else:
                assert pos_name in Mingeo_XLSX_pos_order[4:]
        
        posindex[pos_name] = act_pos_index

        _pos_name = pos_name
        if _pos_name == Mingeo_XLSX_pos_order[-1]: _pos_name = Mingeo_XLSX_pos_order[-2] # For "SC" Scale Check, we use the scale check basic position
        assert _pos_name[0] in dim_utils.TELESCOPE_CARDINAL_DIRECTIONS, pos_name  # East/West/Norht/South
        assert _pos_name[1] in "DU", _pos_name  # Sensor Down / Sensor Up
        
        mag_pos = mag_meas_positions[act_pos_index]
        assert mag_pos.base is mag_meas_positions
        mag_pos["pos_name"] = pos_name                
        mag_pos[dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION] = _pos_name[0]
        mag_pos[dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION]  = DIM_instr.SUSD_2_FLFR( "S"+_pos_name[1] )[1]
        
        if input_file_type == DI_spreadsheet_file_type_enum.Mingeo_manual_XLSX:
            rotated_angle, fixed_angle = manual_xlsx_file__parse_angles ( NM_measpart_act_row_index, series_part )
        elif input_file_type == DI_spreadsheet_file_type_enum.Mingeo_DS_XLSX:
            rotated_angle, fixed_angle = DS_xlsx_file__parse_NMpart_angles(NM_measpart_act_row_index, series_part)
        else:
            rotated_angle, fixed_angle = gen_parse_NMpart_angles(series_part, pos_name, NM_measpart_act_row_index, NM_measpart_rotated_axis_col_index, D_from_decl_series = D_from_file) 
            
        default_datetime = read_in_next_meas_part.last_measpart_datetime + datetime.timedelta(seconds=30) if estiamte_timelabel_if_missing  else np.nan
        measpart_datetime = cell_timestamp_to_datetime( [NM_measpart_act_row_index, timestamp_input_col_index], default_datetime )
      
        if not isinstance(measpart_datetime, datetime.datetime):
            assert False
        
        read_in_next_meas_part.last_measpart_datetime = measpart_datetime
        mag_pos["hc"] = rotated_angle if series_part == "Decl" else fixed_angle
        mag_pos["vc"] = fixed_angle   if series_part == "Decl" else rotated_angle
        mag_pos["timelabel_datetime"] = measpart_datetime 
        mag_pos["timelabel"] = mu.datetime_2_exceltime( mag_pos["timelabel_datetime"]  )
        if series_part == "Decl":
            D_series_timelabels += [ mag_pos["timelabel"] ]
            D_series_indeces += [ act_pos_index ]
        elif series_part == "Incl":
            I_series_timelabels += [ mag_pos["timelabel"] ]
            I_series_indeces += [ act_pos_index ]
        else:
            assert False
        if M_colindex is None:
            mag_pos["M"] = 0.0
        else:
            M = absmeas_table.loc[NM_measpart_act_row_index,M_colindex]
            if type(M) is str: 
                M = float(M) # E.g. : "    0.2" => 0.2
            else:
                assert isinstance(M, (int,float)), 'M: "{}" type:{}'.format(M, type(M) )            
            assert not np.isnan(M) # Asserted, if M cell is empty. Empty M cells may mean of 0 nT ?!
            mag_pos["M"] = float(M)

        NM_measpart_act_row_index += 1
        act_pos_index += 1
        
    # Uses "amdp.dim_instrument_model" to estimate M value for a gives positions
    def calc_M_for_pos(pos_name, M_offset):
        F = gu.measurement_site_list[meas_pos.station_code]["F"]
        vpa_n = amdp.vpa_tn(mag_meas_positions[posindex[pos_name]]['vc']) 
        hpa_n = amdp.hpa_tn(mag_meas_positions[posindex[pos_name]]['hc'], MIRA_mean_for_SU_pos, azimuth_mark_NED )
        M = amdp.dim_instrument_model(F, vpa_n, hpa_n, I_from_file, D_from_file, e = epsilon, d = delta_declination_part, o = M_offset, s = M_scale)
        return M
        
        
    read_in_next_meas_part.last_measpart_datetime  = None
    
    read_in_next_meas_part("Decl")
    read_in_next_meas_part("Decl")
    read_in_next_meas_part("Decl")
    read_in_next_meas_part("Decl")
    
    max_angle_diff  = math.radians(0.6) # In DS files, and in some special Excel tables, the vc angles may slgithly differ
    assert abs(mag_meas_positions[posindex["ED"]]["vc"] -  mag_meas_positions[posindex["WD"]]["vc"]) < max_angle_diff
    assert abs(mag_meas_positions[posindex["EU"]]["vc"] -  mag_meas_positions[posindex["WU"]]["vc"]) < max_angle_diff 
    if (abs(mag_meas_positions[posindex["ED"]]["vc"]-mag_meas_positions[posindex["EU"]]["vc"]) - np.pi) > max_angle_diff :
        assert False
    
    H = gu.measurement_site_list[meas_pos.station_code]["H"]
    D_self_calc, MIRA_mean_for_SU_pos, M_offset_declination_part, delta_declination_part = dim_utils.NM_declination_part_basic_evaluation(mag_meas_positions[D_series_indeces], posindex, start_MIRA_UP,start_MIRA_DN, end_MIRA_UP, end_MIRA_DN,  azimuth_mark_NED, H)

    
    if D_from_file_deg is None:
        assert D_from_file is None
        D_from_file, D_from_file_deg  = D_self_calc, math.degrees(D_self_calc)
    if D_from_file_deg>180: 
        D_from_file_deg = D_from_file_deg - 360.0
    if consistency_checks:
        assert abs(D_from_file_deg - math.degrees(D_self_calc)) < 0.1, "Too big difference between file and calculated declinations: {:.2f} deg / {:.2f} deg".format(D_from_file_deg, math.degrees(D_self_calc))


    # ----------------------------------------------
    # -- Reading in inclination series points  -----
    # ----------------------------------------------

    NM_measpart_act_row_index = NM_I_measpart_first_row_index
    read_in_next_meas_part("Incl")
    read_in_next_meas_part("Incl")
    read_in_next_meas_part("Incl")
    read_in_next_meas_part("Incl")
    
    
    has_scale_check = NM_SC_measpart_row_index is not None
    if has_scale_check:
        NM_measpart_act_row_index = NM_SC_measpart_row_index
        read_in_next_meas_part("Incl", True)
        I_series_timelabels.pop() # Removing "SC" from I series
        I_series_indeces.pop()
        has_scale_check = not np.isnan(mag_meas_positions[posindex["SC"]]['M'])
    else:
        posindex["SC"] = 8
        mag_meas_positions[posindex["SC"]]['M'] = np.nan
        mag_meas_positions[posindex["SC"]]['hc'] = np.nan
        mag_meas_positions[posindex["SC"]]['vc'] = np.nan
        
    assert posindex["SC"] == len(mag_meas_positions) - 1        

    assert abs(mag_meas_positions[posindex["SD"]]["hc"] - mag_meas_positions[posindex["NU"]]["hc"]) < max_angle_diff
    assert abs(mag_meas_positions[posindex["ND"]]["hc"] - mag_meas_positions[posindex["SU"]]["hc"]) < max_angle_diff
    if has_scale_check: assert abs(mag_meas_positions[posindex["SU"]]["hc"] - mag_meas_positions[posindex["SC"]]["hc"]) < max_angle_diff
    if abs( abs(mag_meas_positions[posindex["ND"]]["hc"] - mag_meas_positions[posindex["SD"]]["hc"]) - np.pi )>max_angle_diff:
        assert False
    # assert abs(abs( angle_utils.angle_distance(ND_hc,SD_hc) ) - np.pi) <0.01  ,  "ND_hc:{:.1f}°   SD_hc:{:.1f}°".format( math.degrees(ND_hc), math.degrees(SD_hc) )

    # Safety check: inclination series should start after declination series
    assert mag_meas_positions[posindex["NU"]]["timelabel"] >= mag_meas_positions[posindex["WU"]]["timelabel"], "{} {}".format(mag_meas_positions[posindex["NU"]]["timelabel"],mag_meas_positions[posindex["WU"]]["timelabel"])
    assert mag_meas_positions[posindex["SD"]]["timelabel"] >= mag_meas_positions[posindex["WU"]]["timelabel"]
    assert mag_meas_positions[posindex["ND"]]["timelabel"] >= mag_meas_positions[posindex["WU"]]["timelabel"]
    assert mag_meas_positions[posindex["SU"]]["timelabel"] >= mag_meas_positions[posindex["WU"]]["timelabel"]

    I_self_calc = amdp.Inclination(mag_meas_positions[posindex["NU"]]['vc'], mag_meas_positions[posindex["SD"]]['vc'], mag_meas_positions[posindex["ND"]]['vc'], mag_meas_positions[posindex["SU"]]['vc'])    
    assert abs( I_self_calc - I_from_file) < math.radians(0.5)
    
    epsilon_D  = amdp.Epsilon_D(mag_meas_positions[posindex["ED"]]['hc'], mag_meas_positions[posindex["WD"]]['hc'], mag_meas_positions[posindex["EU"]]['hc'], mag_meas_positions[posindex["WU"]]['hc'], I_self_calc)
    epsilon_I  = amdp.Epsilon_I(mag_meas_positions[posindex["NU"]]['vc'], mag_meas_positions[posindex["SD"]]['vc'], mag_meas_positions[posindex["ND"]]['vc'], mag_meas_positions[posindex["SU"]]['vc'])
    epsilon    = np.mean( [ epsilon_D, epsilon_I]  )
    
   
    if has_scale_check: 
        assert mag_meas_positions[posindex["SC"]]["timelabel"] >= mag_meas_positions[posindex["WU"]]["timelabel"]
        can_calc_SC = True
    else:
        # In residual method, we can miss the scale check measurement, if M values are enough big!
        can_calc_SC = np.sum(np.abs(mag_meas_positions['M'])) > 10.0
        if can_calc_SC: 
            mag_meas_positions = mag_meas_positions[:-1]
            assert len(mag_meas_positions) == 8
            posindex["SC"] = -1
        else:
            # Adding an artificial SC record (with M_scale = 1) as it can be required for post-processing
            mag_meas_positions[posindex["SC"]]['pos_name'] = "SC"
            mag_meas_positions[posindex["SC"]]['hc'] = mag_meas_positions[posindex["SU"]]['hc']
            mag_meas_positions[posindex["SC"]]['vc'] = mag_meas_positions[posindex["SU"]]['vc'] + np.radians(0.2) 
            M_SC = calc_M_for_pos("SC",M_offset_declination_part)
            assert abs(M_SC < 200), f"${M_SC=}"
            mag_meas_positions[posindex["SC"]]['M'] = M_SC
            mag_meas_positions[posindex["SC"]]['timelabel']  = max(mag_meas_positions['timelabel']) + 30.0 / 86400.0
            mag_meas_positions[posindex["SC"]]["timelabel_datetime"] =  mu.exceltime_2_datetime(mag_meas_positions[posindex["SC"]]['timelabel'])

    warning_printed = False
    for mag_meas_position in mag_meas_positions:
        pos_name = mag_meas_position["pos_name"].decode(encoding="utf-8")
        M_calc = calc_M_for_pos(pos_name, M_offset_declination_part)
        M_diff = abs(M_calc - mag_meas_position["M"])
        # We allow negative scale factor in scale check pos:
        if (pos_name == "SC") and abs(M_calc - -mag_meas_position["M"])<M_diff: M_diff = abs(M_calc - -mag_meas_position["M"]) 
        if (M_diff>40) and (not warning_printed):
            warning_printed = True
            cu.warning("\nBig difference between measured and calculated M in pos {:}: measured:{:.1f} nT - calculated:{:.1f} nT".format(pos_name, mag_meas_position["M"], M_calc))
        

    # meas_point_datetimes=(ED_datetime, WD_datetime, EU_datetime, WU_datetime, NU_datetime, SD_datetime, ND_datetime, SU_datetime, SC_datetime)
    # meas_point_timelabels=[]
    # for i in range(len(meas_point_datetimes)): meas_point_timelabels += [ mu.datetime_2_exceltime( meas_point_datetimes[i] ) ]

    # for i in range(len(mag_meas_positions)-1):
    #    assert mag_meas_positions[i]["timelabel"] < mag_meas_positions[i+1]["timelabel"], "Nem hagyományosan folytonos sorrend kezelése nincs még megoldva"
    

    # In same cases, the measurement sequence is difference from the standard one
    mag_part_start_timelabel =  min(mag_meas_positions["timelabel"][:])
    mag_part_end_timelabel   =  max(mag_meas_positions["timelabel"][:])
    
    filter_out_reason = None

    if (datetime_filter is not None):
        if (mag_part_start_timelabel <= mu.datetime_2_exceltime(datetime_filter[0])) or (mag_part_end_timelabel >= mu.datetime_2_exceltime(datetime_filter[1])):
            filter_out_reason = "out of time window"

    # Checking if measurement total time is reasonable
    mag_part_total_time_seconds = (mag_part_end_timelabel - mag_part_start_timelabel)*86400.0
    if  mag_part_total_time_seconds > 3600:
       assert False, f"Measurement took too long:  {mag_part_total_time_seconds} sec"
    if  mag_part_total_time_seconds < 60:
       assert False, f"Measurement took too short time:  {mag_part_total_time_seconds} sec"


    assert len(D_series_timelabels) == len(I_series_timelabels) == 4

    return absmeas_table, filter_out_reason, [mag_part_start_timelabel, mag_part_end_timelabel, I_series_timelabels, D_series_timelabels ],  pillar, [theo_SN, electronics_SN], I_from_file_deg, D_from_file_deg, [start_MIRA_UP, start_MIRA_DN, end_MIRA_UP, end_MIRA_DN], mag_meas_positions, [refmark_descr, azimuth_mark_NED_deg]


# Within some spreadsheet files, there can be more worksheets with different measurements.
# This function lists with measurement data. Listing is based on sheet name mask.
def get_worksheet_names_with_absmeas_data(input_file_type, xlsx_file_name, datetime):
    xls_file = pd.ExcelFile(xlsx_file_name)
    sheet_names = xls_file.sheet_names
    mask = input_file_defs[input_file_type]["worksheet_mask"]
    if (not '?' in mask) and (not '*' in mask):
        # We expect one sheet, with the given name
        assert len(sheet_names) == 1
        assert sheet_names[0] == mask
        return sheet_names
    elif mask == "*":
        # Using all sheets
        return sheet_names
    else:
        mask = datetime.strftime(mask)
        result = []
        for sheet_name in sheet_names:
            if fnmatch.fnmatch(sheet_name, mask):
                result += [sheet_name]
        if len(result)==0:
            assert False, f"{xlsx_file_name=} {input_file_type=} {sheet_names=} {mask=}"
        return result
    

# Reading in Absolute measurement data from XLSX file - one measurement from one file (or from one tabsheet)
def read_absmeas_XLS_files(meas_pos, datelimit_from, datelimit_to, base_path = None, filename_mask = None, F_dat_reader = None, DIM_instr = None, theo_SN_filter = None, input_file_types = supported_input_file_types, consistency_checks = True):

    # assert (DIM_instr is None) or isinstance(DIM_instr, globals()['DIM_instr'].ADS_model_class)
    import dim_measurement_loader
    from dim_measurement_loader import COMPNAME__PROGRAM_MAINTYPE_INDEX, COMPNAME__TIMESTAMP, COMPNAME__POS_ID, COMPNAME__MAG, COMPNAME__MEASURED_HC_RAD, COMPNAME__MEASURED_VC_RAD, COMPNAME__FILE_ABMEAS_WIN_INDEX, COMPNAME__INPUT_FILE_INDEX, DATAWINDOW_COMPNAME__M_RECORD_INDECES, DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES, COMPNAME__AZIMARK_INDEX,  COMPNAME__PILLAR_INDEX, COMPNAME__AZIMARK_WGS84, COMPNAME__TELESCOPE_CARDINAL_DIRECTION

    def unexpected_cell_data_errmsg(absmeas_table, cell_def, header_cell_relative_index = None):
        cell_data = absmeas_table.loc[ xlsx_cells_def[cell_def][0], xlsx_cells_def[cell_def][1] ]
        if header_cell_relative_index is not None:
            header_cell_data =  '"'+absmeas_table.loc[ xlsx_cells_def[cell_def][0], xlsx_cells_def[cell_def][1]+header_cell_relative_index ]+'"'
            header_cell_data += ":"
        else:
            header_cell_data = ""

        return "Unexpected cell data: {}\"{}\"".format(header_cell_data,cell_data)
        assert False
        sys.exit()


    # recalc_I = []
    # recalc_D = []
    #for i in recalc_methods:
    #    recalc_I += [[]]
    #    recalc_D += [[]]

    # read_all_pillars_at_once =  or True

    # if read_all_pillars_at_once:
    #     XLS_path  = absmeas_xlsx_path
    #     pillars_in_path = pillars_in_absmeas_xlsx_path
    #     loop_count = 1
    #     # act_pillars = pillars
    # else:
    #     XLS_path  = None
    #     act_pillar = None
    #     loop_count = len(pillars)


    file_params_list = []
    file_start_datetime_list = []
    
    # All measurements from XLSX files is supposed to be in this type
    program_maintype_str_list = ["MAN-RESID"]    
    
    use_default_filename_mask = filename_mask is None
    if not use_default_filename_mask: 
        # amdp run gets here
        pillars_subdir = None
        if isinstance(filename_mask , (list,tuple)):
            assert len(filename_mask) == 1
            filename_mask = filename_mask[0]
        
    pillars = meas_pos.pillars
    print("---- STEP-1 Creating DIM measurement file list, sorted by date -----")
    nof_abs_measurements = 0
    for input_file_type in input_file_types:
        
        input_file_def = input_file_defs[input_file_type]
        
        if use_default_filename_mask:
           act_filename_mask  = input_file_def['default-filename-mask']  # System default mask 
        else:
           act_filename_mask = base_path + filename_mask  # User supplied mask       
        
        if (len(pillars)>1) and ( PILLAR_IN_FILENAME_MASK in act_filename_mask ) :
            pillar_search_loop_count = len(pillars)
        else:
            pillar_search_loop_count = 1
    
        for pillar_search_index in range(pillar_search_loop_count):
            
            if use_default_filename_mask: pillar_subdir = input_file_def['pillar-subdir'] 
            
            if pillar_search_loop_count>1:
                # If different pillars are in different subdirectories
                act_pillar = pillars[pillar_search_index]
                if pillars_subdir is not None: pillars_subdir  = pillars_subdir.replace(PILLAR_IN_FILENAME_MASK, act_pillar) 
                pillars_in_path = [act_pillar]
            else:
                pillars_in_path = pillars
                act_pillar = pillars[0]
                pass
        
                filename_mask  = filename_mask.replace(STATION_IN_FILENAME_MASK,meas_pos.station_code)
                filename_mask  = filename_mask.replace(INSTR_IN_FILENAME_MASK,"*")
    
                #if not read_all_pillars_at_once:
                filename_mask = filename_mask.replace(PILLAR_IN_FILENAME_MASK, act_pillar)
    
            # filename_datetime_part_len = len(datetime.datetime.now().strftime(filename_mask))
            
            # Date in file name mask:
            if not "%Y" in filename_mask:
                assert not "%m" in filename_mask
                assert not "%d" in filename_mask
                filename_date_len = None
                act_measurement_date_is_read_from_cell = True
            else:
                filename_date_len = filename_mask.index("%d") - filename_mask.index("%Y") + 4
                filename_datetime_mask = filename_mask[ filename_mask.index("%Y"): filename_mask.index("%d")+2]  # Date parsing part            
                filename_mask = filename_mask.replace("%Y","????").replace("%m","??").replace("%d","??")
                act_measurement_date_is_read_from_cell = def_measurement_date_is_read_from_cell
            
            full_path = base_path
            if pillars_subdir is not None: full_path = os.path.join(full_path, pillars_subdir)
            full_path = full_path.replace(STATION_IN_FILENAME_MASK, meas_pos.station_code)
            print('Reading files "{}" from dir "{}"...'.format(filename_mask,full_path), end='')
            # XLSX_filelist_orig = glob.glob( os.path.join(full_path, filename_mask ) )
            # XLSX_filelist_orig = sorted(XLSX_filelist) # Sorting file name to have increasing dates:
            
            # Testing a new method:
            XLSX_filelist, used_masks = mu.filenamemask_to_filenames(full_path, filename_mask)
            XLSX_filelist = sorted(XLSX_filelist)
            # assert len(XLSX_filelist_orig) == len(XLSX_filelist)
            # for i in range(len(XLSX_filelist)):
            #     assert XLSX_filelist[i] == full_path + XLSX_filelist_2[i], "{} <> {}".format(XLSX_filelist[i], XLSX_filelist_2[i])
            
            nof_excel_data_files=len(XLSX_filelist)
            if nof_excel_data_files==0:
                print ("No data file found!")
                continue
            else:
                 print("  found {:} files.".format(nof_excel_data_files))
    
            
        # Now determining dates for all files
            for i,filename_only in enumerate(XLSX_filelist):
                full_filename = os.path.join(full_path, filename_only)
                
                print("{:}/{:} ({:})  \"{:}\": ".format(i+1,nof_excel_data_files,len(file_params_list)+1, filename_only), end='')
                if filename_date_len:
                    date_start_pos = filename_only.index("_20")+1  # There are two file name formats: "MNY_2025-..." and "MNY_A1_2025-..."
                    date_str = filename_only[date_start_pos:date_start_pos+filename_date_len]
                    # date_str = filename_only[:filename_datetime_part_len],filename_datetime_maks)
                    date_from_filename=datetime.datetime.strptime(date_str,filename_datetime_mask)
                    assert (date_from_filename.hour == 0) and (date_from_filename.minute == 0)
                else:
                    # File name does not have date
                    date_from_filename = None
    
                if not act_measurement_date_is_read_from_cell and (date_from_filename is not None):
                    # 'datelimit_from' and 'datelimit_to' can have hour/minute/second parts, while 'date_from_filename' has only date part
                    if (date_from_filename.date() < datelimit_from.date()) or  (date_from_filename.date() > datelimit_to.date()) :
                        print("out of time window - skipping")
                        continue

                worksheet_names_with_dim_data = get_worksheet_names_with_absmeas_data(input_file_type, full_filename, date_from_filename)
                first_worksheet_name = worksheet_names_with_dim_data[0]
                
                absmeas_table, filter_out_reason, [start_timelabel, end_timelabel, I_series_timelabels, D_series_timelabels], pillar, [Theo_SN, electronics_SN], I_deg_not_used, D_deg_not_used,  [start_MIRA_UP, start_MIRA_DN, end_MIRA_UP, end_MIRA_DN], mag_meas_positions, [refmark_descr, azimuth_mark_NED_deg]  = read_absmeas_XLS_file( input_file_type, full_filename, first_worksheet_name, meas_pos, date_from_filename, theo_SN_filter, [datelimit_from,datelimit_to], DIM_instr = DIM_instr, consistency_checks = consistency_checks) # , instr, theo_SN_list
                del absmeas_table
    
                # include_all_pillars = pillars[0] == '*'
                if (filter_out_reason is None) and (not pillar in pillars):  # and  (not include_all_pillars)  
                    filter_out_reason = "excluded pillar \"{}\"".format(pillar)
    
                if filter_out_reason is not None:
                     print(filter_out_reason, "-> not using")
                     continue
    
                if pillars_in_path is not None:
                    assert ( pillar in pillars_in_path), f"Expected pillar: {pillars_in_path} / Pillar from file:{pillar}"
                    
                # Determining instrument name, that is plotted on baseline plot 
                if input_file_def == DI_spreadsheet_file_type_enum.Mingeo_manual_XLSX:
                    instr_name = filename_only.split("_")[-1]  # File name fields are separated by '_', and the last part is the instrument, eg. "REFE" (not always!!!)
                else:
                    instr_name = Theo_SN
                if instr_name == "354459": instr_name= "REF"
                    
    
                file_params_list += [ {"file_type":input_file_type, "full_file_name":full_filename, "file-name":filename_only, "start-date":start_timelabel, "end-date":end_timelabel, "instr-name":instr_name, "pillar-index": pillars.index(pillar), "date_from_filename":date_from_filename }  ] 
                
                file_start_datetime_list += [start_timelabel]
                
                nof_abs_measurements += len(worksheet_names_with_dim_data)
    
                print("within time window -> using")

    if len(file_params_list) == 0:
        return "No DIM measurement file found"

    print("---- STEP-2 Sorting files by date -----")

    sort_order = np.argsort( np.array(file_start_datetime_list) )
    # file_params_list = np.array(file_params_list)[sort]
    # file_params_list = list(file_params_list)
    file_params_list  = [  file_params_list[i] for i in sort_order ]
    del file_start_datetime_list
    
    for i in range(len(file_params_list)-1):
        if file_params_list[i]["start-date"]  >= file_params_list[i+1]["start-date"]:
            fname1=file_params_list[i]["file-name"]
            fname2=file_params_list[i+1]["file-name"]
            mu.warning("Measurement time jump or repeated measurement: {}-{}".format(fname1,fname2 ) )
        if file_params_list[i]["end-date"]  >= file_params_list[i+1]["start-date"]:
            mu.warning("Time overlap: {}  - {}".format(file_params_list[i]["file-name"], file_params_list[i+1]["file-name"]))
    
    print("---- STEP-3 Processing DIM measurement files -----")

    file_absmeas_window_data = np.empty( nof_abs_measurements, dtype = 
            dim_measurement_loader.FILE_ABSMEAS_WINDOW_DATA_DTYPE + [ 
            ("starttime", float), 
            ("endtime" , float), 
            ("reftime", float),            
            ("Excel-I", float ), 
            ("Excel-D", float ), 
            ("theo-SN", "S16"),
            ("electronics-SN", "S16"),  ]  
    #      # ("inclseries_time", "O"), # it will store a Python object reference (pointer) to an ndarray
    #      # ("declseries_time", "O"), 
        
        )
                                        
    refmark_data_       = dict()
    refmark_data_index_ = dict()
    abs_meas_window_index, M_data_index, refmark_data_index_["FL"], refmark_data_index_["FR"] = 0, 0, 0, 0
    
    dtype = [
           (COMPNAME__FILE_ABMEAS_WIN_INDEX, int) ,
           (COMPNAME__MEASURED_HC_RAD, float),
           (COMPNAME__MEASURED_VC_RAD, float),
           (COMPNAME__TIMESTAMP, float),
           (COMPNAME__POS_ID, int),
           (COMPNAME__INPUT_FILE_INDEX, int),
           (dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION, dim_measurement_loader.COMPTYPE__TELESCOPE_CARDINAL_DIRECTION),
           (dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION,  dim_measurement_loader.COMPTYPE__TELESCOPE_FACE_POSITION),
        ] 
    M_data =       np.empty( nof_abs_measurements * NR_OF_DIM_MEAS_MAIN_POSITIONS, dtype = dtype + [ (dim_measurement_loader.COMPNAME__MAG,  float)] )
    for FLFR in ["FL", "FR"]:
        refmark_data_[FLFR] = np.empty( nof_abs_measurements * NR_OF_DIM_REFMARK_MEAS_POSITIONS // 2, dtype = dtype + [  (COMPNAME__AZIMARK_INDEX, int), (COMPNAME__PILLAR_INDEX, int), (COMPNAME__AZIMARK_WGS84, float) ]  )   

    refmark_list = []
    for input_file_index, file_params in enumerate(file_params_list):
        pillar_index = file_params["pillar-index"]
        full_filename = file_params["full_file_name"]
        file_type = file_params["file_type"] 
        instr_name = file_params["instr-name"]
        filename_only = file_params["file-name"]
        file_datetime = file_params["date_from_filename"] 
        worksheet_names = get_worksheet_names_with_absmeas_data(file_type, full_filename, file_datetime)
        for worksheet_name in worksheet_names:
        
            absmeas_table, filter_out_reason, [start_timelabel, end_timelabel, I_series_timelabels, D_series_timelabels], pillar, [Theo_SN, electronics_SN], I_from_file_deg, D_from_file_deg,  [start_MIRA_UP, start_MIRA_DN, end_MIRA_UP, end_MIRA_DN], mag_meas_positions, [refmark_descr, azimuth_mark_NED_deg ] = read_absmeas_XLS_file(file_type, full_filename, worksheet_name, date_from_filename = file_datetime, meas_pos = meas_pos, theo_SN_filter = None , datetime_filter = None, DIM_instr = DIM_instr, consistency_checks = consistency_checks )
            assert filter_out_reason is None
    
            ref_timelabel   =  np.mean([start_timelabel,end_timelabel])
            

            #if len(file_absmeas_window_data)<nof_abs_measurements:
            #    file_absmeas_window_data = np.resize( file_absmeas_window_data, nof_abs_measurements + 10 ) 
            assert len(file_absmeas_window_data)>=abs_meas_window_index + 1
                
            act_abmeas_window_data = file_absmeas_window_data[abs_meas_window_index]
            
            act_abmeas_window_data[COMPNAME__PROGRAM_MAINTYPE_INDEX] = 0  # Index within 'program_maintype_str_list' 
            act_abmeas_window_data["file-index"] =  input_file_index
            act_abmeas_window_data["Excel-I"] = math.radians(I_from_file_deg)
            act_abmeas_window_data["Excel-D"] = math.radians(D_from_file_deg)
            
            act_abmeas_window_data["starttime"] = start_timelabel
            act_abmeas_window_data["endtime"]   = end_timelabel
            act_abmeas_window_data["reftime"]   = ref_timelabel            
            
            # act_abmeas_window_data["inclseries_time"] = np.array(I_series_timelabels, dtpye = "f8")
            # act_abmeas_window_data["declseries_time"] = np.array(D_series_timelabels, dtpye = "f8")
            
            act_abmeas_window_data["pillar-index"] = pillar_index
            
            act_abmeas_window_data["theo-SN"] =   Theo_SN 
            act_abmeas_window_data["electronics-SN"] = electronics_SN 
    
            print("{:}/{:} \"{:}\": DateTime={:} I={:.4f}°   D={:.4f}°".format(input_file_index, len(file_params_list), filename_only, mu.exceltime_2_datetime(start_timelabel), I_from_file_deg, D_from_file_deg))
    
            #for i_recalc_method,recalc_method in enumerate(recalc_methods):
            #    _recalc_I, _recalc_D = post_reprocess_manual_absmeas(recalc_method)
            #    recalc_I[i_recalc_method] += [_recalc_I]
            #    recalc_D[i_recalc_method] += [_recalc_D]
            
            act_abmeas_window_data[COMPNAME__FILE_ABMEAS_WIN_INDEX] = abs_meas_window_index
            # act_abmeas_window_data["M_start_record_index"] = M_data_index 
            # act_abmeas_window_data["M_end_record_index"]   = M_data_index+NR_OF_DIM_MEAS_MAIN_POSITIONS-1
            # act_abmeas_window_data["variation_idx"] = 0
            # act_abmeas_window_data["name"] = 


            # Preparing returned arrays:            
            if meas_pos.station_code == "MNY":
                assert pillar in [ "A1", "A2" ]
                assert pillar_index == meas_pos.pillars.index(pillar), f"{pillar} != {meas_pos.pillars}[{pillar_index}]"
                refmark_name = "TowerLeft" if (refmark_descr[:11] == "Water Tower") else refmark_descr
                # assert refmark_descr[:11] == "Water Tower" , unexpected_cell_data_errmsg(absmeas_table,  'refmark_descr', -2)
            else:
                refmark_name = "X:{}".format(azimuth_mark_NED_deg)
                
            if not refmark_name in refmark_list: refmark_list += [refmark_name]
            MIRA_UP_horiz_angle = DIM_instr.get_horizontal_vc_angle("SU","rad")
            MIRA_DN_horiz_angle = DIM_instr.get_horizontal_vc_angle("SD","rad")
            time_ordered_mag_meas_positions = mag_meas_positions[ np.argsort(mag_meas_positions, order="timelabel") ]
            if F_dat_reader is not None:
                # Reading F from a data file:
                timelabel_sink,F_variation, sink = F_dat_reader.get_data("RAW_INTERPOLATED", time_ordered_mag_meas_positions['timelabel'])

            M_reading_record_indeces = list( range(M_data_index, M_data_index + NR_OF_DIM_MEAS_MAIN_POSITIONS ) )            
            M_data[COMPNAME__INPUT_FILE_INDEX]       [M_data_index:M_data_index+NR_OF_DIM_MEAS_MAIN_POSITIONS] =  input_file_index
            M_data[COMPNAME__FILE_ABMEAS_WIN_INDEX]  [M_data_index:M_data_index+NR_OF_DIM_MEAS_MAIN_POSITIONS] =  abs_meas_window_index
            # raw_data[COMPNAME__PROGRAM_MAINTYPE_INDEX] [raw_data_index:raw_data_index+RAWDATA_RECORDS_PER_ABSMEAS] = 0  # Index within 'program_maintype_str_list'

            REFPOS_MEASUREMENT_RECORDS_POS_ID = -1
            # refmark_reading_record_indeces = list( range(refmark_data_index_FLFR,refmark_data_index_FLFR+NR_OF_DIM_REFMARK_MEAS_POSITIONS)  )                
            assert DIM_instr.SUSD_2_FLFR("SU") == "FL" # The lines below suppose the standard sensor position
            refmark_reading_record_indeces = {"FL": [refmark_data_index_["FL"]+0, refmark_data_index_["FL"]+1 ], "FR": [refmark_data_index_["FR"]+0, refmark_data_index_["FR"]+1] }

            for FLFR in "FL", "FR":
                refmark_data_[FLFR] [COMPNAME__INPUT_FILE_INDEX][refmark_reading_record_indeces[FLFR]] = input_file_index
                refmark_data_[FLFR] [COMPNAME__FILE_ABMEAS_WIN_INDEX][refmark_reading_record_indeces[FLFR]] = abs_meas_window_index
                refmark_data_[FLFR] [COMPNAME__POS_ID][refmark_reading_record_indeces[FLFR]] = REFPOS_MEASUREMENT_RECORDS_POS_ID
                refmark_data_[FLFR] [COMPNAME__AZIMARK_INDEX][refmark_reading_record_indeces[FLFR]] = refmark_list.index(refmark_name)
                refmark_data_[FLFR] [COMPNAME__PILLAR_INDEX][refmark_reading_record_indeces[FLFR]] = pillar_index  # Note: this data is reduntant, as it can be read out from file data using 'act_abmeas_window_data["file-index"]'
                refmark_data_[FLFR] [COMPNAME__AZIMARK_WGS84][refmark_reading_record_indeces[FLFR]] = np.radians(azimuth_mark_NED_deg) # This is redundant also in amdp
                refmark_data_[FLFR] [COMPNAME__TELESCOPE_CARDINAL_DIRECTION] = dim_utils.UNKNOWN_CARDINAL_DIRECTION
            
            refmark_data_["FL"][COMPNAME__TIMESTAMP][refmark_data_index_["FL"]+0] = time_ordered_mag_meas_positions[0]["timelabel"] - 2 * REFMARK_READING_ESTIMATED_DURATION__SEC / 86400.0
            refmark_data_["FL"][COMPNAME__MEASURED_HC_RAD][refmark_data_index_["FL"]+0] = start_MIRA_UP
            refmark_data_["FL"][COMPNAME__MEASURED_VC_RAD][refmark_data_index_["FL"]+0] = MIRA_UP_horiz_angle
            
            refmark_data_["FR"][COMPNAME__TIMESTAMP][refmark_data_index_["FR"]+0] = time_ordered_mag_meas_positions[0]["timelabel"] - 1 * REFMARK_READING_ESTIMATED_DURATION__SEC / 86400.0
            refmark_data_["FR"][COMPNAME__MEASURED_HC_RAD][refmark_data_index_["FR"]+0] = start_MIRA_DN
            refmark_data_["FR"][COMPNAME__MEASURED_VC_RAD][refmark_data_index_["FR"]+0] = MIRA_DN_horiz_angle
            
            refmark_data_["FL"][COMPNAME__TIMESTAMP][refmark_data_index_["FL"]+1] = time_ordered_mag_meas_positions[0]["timelabel"] + 1 * REFMARK_READING_ESTIMATED_DURATION__SEC / 86400.0
            refmark_data_["FL"][COMPNAME__MEASURED_HC_RAD][refmark_data_index_["FL"]+1] = end_MIRA_UP
            refmark_data_["FL"][COMPNAME__MEASURED_VC_RAD][refmark_data_index_["FL"]+1] = MIRA_UP_horiz_angle
            
            refmark_data_["FR"][COMPNAME__TIMESTAMP][refmark_data_index_["FR"]+1] = time_ordered_mag_meas_positions[0]["timelabel"] + 2 * REFMARK_READING_ESTIMATED_DURATION__SEC / 86400.0
            refmark_data_["FR"][COMPNAME__MEASURED_HC_RAD][refmark_data_index_["FR"]+1] = end_MIRA_DN
            refmark_data_["FR"][COMPNAME__MEASURED_VC_RAD][refmark_data_index_["FR"]+1] = MIRA_DN_horiz_angle
            
            
            refmark_data_index_["FL"] += NR_OF_DIM_REFMARK_MEAS_POSITIONS // 2
            refmark_data_index_["FR"] += NR_OF_DIM_REFMARK_MEAS_POSITIONS // 2
            
            assert len(time_ordered_mag_meas_positions) == NR_OF_DIM_MEAS_MAIN_POSITIONS
            M_data[COMPNAME__POS_ID][M_reading_record_indeces] = np.arange(0, NR_OF_DIM_MEAS_MAIN_POSITIONS)
            # M_data[COMPNAME__IS_M_RECORD][M_reading_record_indeces] = True
            # M_data[COMPNAME__AZIMARK_INDEX][M_reading_record_indeces] = INV_REFMARK_INDEX
            M_data[COMPNAME__TIMESTAMP][M_reading_record_indeces] = time_ordered_mag_meas_positions["timelabel"]
            M_data[COMPNAME__MAG][M_reading_record_indeces] = time_ordered_mag_meas_positions["M"]
            M_data[COMPNAME__MEASURED_HC_RAD][M_reading_record_indeces] = time_ordered_mag_meas_positions["hc"]
            M_data[COMPNAME__MEASURED_VC_RAD][M_reading_record_indeces] = time_ordered_mag_meas_positions["vc"]
            M_data[dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION][M_reading_record_indeces] = time_ordered_mag_meas_positions[dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION]
            M_data[dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION][M_reading_record_indeces] = time_ordered_mag_meas_positions[dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION]
            
            M_data_index += NR_OF_DIM_MEAS_MAIN_POSITIONS
            
            #for pos_id, mag_meas_pos  in enumerate(time_ordered_mag_meas_positions):
            #    assert not np.isnan(mag_meas_pos["M"])
            #    raw_data[raw_data_index][COMPNAME__TIMESTAMP] = mag_meas_pos["timelabel"] 
            #    raw_data[raw_data_index][COMPNAME__MAG] =  mag_meas_pos["M"]
            #    raw_data[raw_data_index][COMPNAME__MEASURED_HC_RAD] = mag_meas_pos["hc"]
            #    raw_data[raw_data_index][COMPNAME__MEASURED_VC_RAD] = mag_meas_pos["vc"]
            #    raw_data[raw_data_index][COMPNAME__AZIMARK_INDEX] = INV_REFMARK_INDEX
            #    raw_data_index += 1
            #    # raw_data[raw_data_index]["F"] =  None if F_dat_reader is None else F_variation[pos_id,0] 
                
                
            
            act_abmeas_window_data[DATAWINDOW_COMPNAME__M_RECORD_INDECES]  = M_reading_record_indeces
            
            act_abmeas_window_data[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FL"]] =  refmark_reading_record_indeces["FL"]
            act_abmeas_window_data[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FR"]] =  refmark_reading_record_indeces["FR"]
                
            abs_meas_window_index += 1                            
            del absmeas_table

    assert abs_meas_window_index == nof_abs_measurements
    # file_absmeas_window_data = np.resize( file_absmeas_window_data, nof_abs_measurements )         
    if M_data is not None:
        assert M_data_index == len(M_data)
        assert refmark_data_index_["FL"] == len(refmark_data_["FL"])
        assert refmark_data_index_["FR"] == len(refmark_data_["FR"])

    print("Total nr. of absolute measurements within the period: {}".format( nof_abs_measurements ))
    
    # Printing all read abs. meas in a table
    b_more_pillars = len(pillars)>1
    table_header = ["Date"]
    if b_more_pillars: table_header += ["Pillar\nindex"]
    table_header += [ "SN", "Duration(min)" ]
    table = [   ]
    for rec in file_absmeas_window_data:
        # mid_date_str = "{}".format( mu.exceltime_2_datetime(result["ref-time"][i]) )
        start_date_str = "{}".format( mu.exceltime_2_datetime(rec["starttime"]) )
        duration_min = (rec["endtime"] - rec["starttime"]) * 1440.0
        record = [ start_date_str ]
        if b_more_pillars: record +=  [ pillars[rec["pillar-index"]] ]
        record += [ "{} / {}".format(rec["theo-SN"].decode("UTF-8"), rec["electronics-SN"].decode("UTF-8")),  "{:.1f}".format(duration_min) ]
        table += [ record ]

    print(tabulate.tabulate(table, headers=table_header, tablefmt="pretty"))

    return (file_params_list, file_absmeas_window_data, M_data, refmark_data_, refmark_list, program_maintype_str_list)


