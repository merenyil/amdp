# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 19:43:33 2026

@author: LaszloM

DIM measurement files' XLSX to MagPy converter

"""

import sys

sys.path.insert(0, "../") # In amdp project, this file is in "utilities" directory, while
                       # 'absmeas_xslx_import' is in the main "amdp" directory
                       

import absmeas_xlsx_import
import geomag_utils as gu
import datetime
import dim_instr
import cli_utils as cu
import misc_utils as mu
import numpy as np
import math
import os

from dim_measurement_loader import COMPNAME__PROGRAM_MAINTYPE_INDEX, COMPNAME__TIMESTAMP, COMPNAME__POS_ID, COMPNAME__MAG, COMPNAME__MEASURED_HC_RAD, COMPNAME__MEASURED_VC_RAD, COMPNAME__FILE_ABMEAS_WIN_INDEX, COMPNAME__INPUT_FILE_INDEX, DATAWINDOW_COMPNAME__M_RECORD_INDECES, DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES, COMPNAME__AZIMARK_INDEX,  COMPNAME__PILLAR_INDEX, COMPNAME__AZIMARK_WGS84, COMPNAME__TELESCOPE_CARDINAL_DIRECTION, COMPNAME__TELESCOPE_FACE_POSITION


piers_data_file = "../piers-data.yaml"

gu.pillar_data_init(piers_data_file)

meas_pos = gu.measurement_pos_class("MNY")
meas_pos.pillars = [ "A1", "A2" ]

datelimit_from = datetime.datetime(2026,10,1)
datelimit_to   = datetime.datetime(2026,10,16)

base_path = r"A:\DAT\MNY\AbsMeas"
filename_mask = "MNY_*%Y-%m-%d_*.xlsx"
xlsx_file_types = [ "Mingeo-manual-XLSX" ]

base_path = r"A:\DAT\MNY\AbsMeas/DS"
xlsx_file_types = [ "Mingeo-DS-XLSX" ]
filename_mask = "*MNY-null_2026-10*xlsx"
M_scale_corr = -1

cu.wait_for_enter_at_warnings = False  # In DS files, "D" or "I" values are many time wrong, resulting in False warnings/errors
consistency_checks = False

output_file_dir = r"A:\DAT\MNY\AbsMeas\MagPy-format"

F_dat_reader = None
DIM_instr = dim_instr.DIM_instrument_model_class( {}, "DIM" )


spreadsheet_results = absmeas_xlsx_import.read_absmeas_XLS_files(meas_pos, datelimit_from, datelimit_to, base_path, filename_mask, F_dat_reader, DIM_instr, theo_SN_filter = None, input_file_types = xlsx_file_types, consistency_checks = consistency_checks)

assert type(spreadsheet_results) is not str, "Error reading spreadsheet files: {}".format(spreadsheet_results)
(input_file_data,  file_absmeas_window_data, M_readings, refmark_readings, azimuth_mark_list, program_maintype_str_list) = spreadsheet_results 

for windows_data in file_absmeas_window_data:
    pillar_index  = windows_data['pillar-index']
    M_indeces     = windows_data[DATAWINDOW_COMPNAME__M_RECORD_INDECES]
    M             = M_readings[M_indeces]
    refmark_indeces_FL,refmark_indeces_FR = windows_data[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FL"]], windows_data[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FR"]]
    theo_SN       = windows_data["theo-SN"].decode('utf-8')  
    start_date    = windows_data["starttime"]
    start_date    = mu.exceltime_2_datetime(start_date)
    
    observer  = "Mingeo personel"
    
    azimuth_mark_dir_Deg = np.degrees( refmark_readings['FL'][ refmark_indeces_FL[0] ] [COMPNAME__AZIMARK_WGS84] )
    
    
    filename = start_date.strftime("%Y%m%d_%H%M%S_{}_{}.txt".format(meas_pos.pillars[pillar_index], meas_pos.station_code) )
    
    with open(os.path.join(output_file_dir, filename) , 'w', encoding='utf-8') as f:
        f.write("# MagPy Absolutes\n")
        f.write("# Abs-Observer:	{}\n".format(observer))
        f.write("# Abs-Theodolite:	{}\n".format(theo_SN))
        f.write("# Abs-TheoUnit:	deg\n")
        f.write("# Abs-FGSensor:	Sensor\n")	
        f.write("# Abs-AzimuthMark:	{:.4f}\n".format(azimuth_mark_dir_Deg))
        f.write("# Abs-Pillar:	{}\n".format(meas_pos.pillars[pillar_index]))
        f.write("# Abs-Scalar:	/\n")
        f.write("# Abs-InputDate: 	{}\n".format(start_date.strftime("%Y-%m-%d")))
        f.write("# Abs-Temperature:	20.0\n")
        f.write("# Abs-Notes: 	    \n")
        
        f.write("Miren\n")
        
        assert len(refmark_indeces_FL) == 2
        assert len(refmark_indeces_FR) == 2
        
        SD_start = math.degrees(refmark_readings['FR'][refmark_indeces_FR[0]][COMPNAME__MEASURED_HC_RAD])
        SU_start = math.degrees(refmark_readings['FL'][refmark_indeces_FL[0]][COMPNAME__MEASURED_HC_RAD])
        SD_end   = math.degrees(refmark_readings['FR'][refmark_indeces_FR[1]][COMPNAME__MEASURED_HC_RAD])
        SU_end   = math.degrees(refmark_readings['FL'][refmark_indeces_FL[1]][COMPNAME__MEASURED_HC_RAD])
        
        
        f.write("{:.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}\n".format(SD_start, SD_start, SU_start, SU_start, SD_end, SD_end, SU_end, SU_end  ) )
    
        f.write("Positions:\n")
        
        pos_order = ( "WU", "EU", "WD", "ED", "NU", "SD", "ND", "SU", "SU_ScaleCheck" )
        pos_order = ( "EU", "WU", "ED", "WD", "NU", "SD", "ND", "SU", "SU_ScaleCheck" )
        
        for pos in pos_order:
            cardinal_dir = pos[0]            
            face_pos = { "U":"L", "D": "R" }[pos[1]]
            M_pos = M[ M[COMPNAME__TELESCOPE_CARDINAL_DIRECTION]  == cardinal_dir ]
            M_pos = M_pos[ M_pos[COMPNAME__TELESCOPE_FACE_POSITION]  == face_pos ] 
            
            repeat_measurement = 2
            if len(M_pos) != 1:
                # This should be the scale check point!!!!
                assert len(M_pos) == 2
                if pos == "SU_ScaleCheck":
                    repeat_measurement = 1
                else:
                    assert (pos == "SU")
                
                
            M_pos = M_pos[0]
            timestamp = mu.exceltime_2_datetime(M_pos[COMPNAME__TIMESTAMP])
            for M_repeat in range(repeat_measurement):
                f.write("{} {:9.4f} {:9.4f}  {:6.2f}\n".format( timestamp.strftime("%Y-%m-%d_%H:%M:%S"), math.degrees(M_pos[COMPNAME__MEASURED_HC_RAD]), math.degrees(M_pos[COMPNAME__MEASURED_VC_RAD]), M_pos[COMPNAME__MAG]  * M_scale_corr) )
        
        f.close()
            
            
            

                 
        
        
