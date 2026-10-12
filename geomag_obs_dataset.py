# -*- coding: utf-8 -*-
"""

geomag_obs_dataset.py

Created on Mon Jul 27 21:54:27 2026

@author: Laszlo Merenyi

"""

import os
import numpy as np

# import geomag_utils_preferences as preferences
# import measurement_sites
import geomag_utils as gu
import misc_utils as mu
import dat_file_reader
import matplotlib.pyplot as plt
import plot_utils
import yaml
# from enum import Enum
import angle_utils
from numpy.lib import recfunctions as rfn
import geomag_utils_preferences as preferences

mag_obs_dataset_def_list= {}
# new_mag_obs_data_reader_id = 0
mag_obs_dataset_defs_yaml_path = None


# Storing the latest 3 magnetic data file reader in 'latest_used_dat_file_readers' to fasten the data file reading.
_latest_used_dat_file_readers=[ [ None, None ], [ None, None ], [ None, None]    ]


# from su:    
#class F_DATA_SOURCES(enum.IntEnum):
#    VARIOMETER_DATA = 0 # "FDATA-FROM-VARDATA"  # Flag of 'F_data_def" parameter of 'dim_mag_measurement_process' function
#    SITES_DEFAULT_F = enum.auto() #    "sites_default_F"
#    DIM_DATA        = enum.auto() # "dim_data"
    
    
# Extension of dat_file_reader_class for DIM data process.
# - reading config from a YAML record
# - Separated mag_comps and aux_comps


class mag_obs_data_file_reader(dat_file_reader.dat_file_reader_class):


     def __init__(self, reader_def, parent):
        # global new_reader_id
        assert type( reader_def ) == dict
        self.parent_mag_obs_data_set = parent
        self.aux_comps = []
        
        
        full_filename_mask = reader_def["path"]
        base_path_mask, filename_mask = os.path.split(full_filename_mask)
        
        file_format = reader_def["file-format"]
        assert file_format in dat_file_reader.SUPPORTED_FILE_FORMATS
        
        interval_def = reader_def["interval"]
        file_data_name = reader_def.get("name", parent.dataset_name)
        
        components_def = reader_def["components"]  # Can be an list of magnetic components (e.g. ["X","Y","Z"], or a full record definition, with timestamp components, and optional skipped components. see dat_file_reader.py
        
        
        if file_format  == "IAGA-2002":
            assert interval_def in ['1min','1sec']  # ToDo: extend the list
            # if self.file_name_format=="*" # File name is automatically determined here:
            # assert data_type in ['q', 'p', 'v', 'd'], f"{data_type=}"  # ToDo: extend the list
            # IAGA-2002 file name format:
            # file_name_mask = IAGA_code.lower()+"{:4d}{:02d}{:02d}"+ data_type + interval_def + "." + interval_def
            #if mag_comps in [ "XYZ", "HDZ" ]:
            #    mag_comps = [mag_comps[0], mag_comps[1], mag_comps[2]]
            #    comps_columns = [ 3, 4, 5 ]
            #elif mag_comps in [ "XYZF", "HDZF" ]:
            #    mag_comps = [mag_comps[0], mag_comps[1], mag_comps[2], "F"]
            #    comps_columns = [ 3, 4, 5, 6]
             
            if len(components_def) == 1  and  (components_def == 'F'):
                comps_columns = [ 6 ] # Reading only "F" from the IAGA-2002 file
            else:
                assert len(components_def) <= 4
            comps_columns = [ 3, 4, 5, 6 ] [0:len(components_def)  ] # Typically 3 or 4 components                
        else:
            comps_columns = None # will be determined automatically
   

        # May be added: FILE_COLINDEX_HOUR = self.hour_colindex 
        super().__init__(base_path_mask, filename_mask, interval_def, components_def, comps_columns, file_format = file_format, data_descr = file_data_name)
        self.mag_comps = self.comp_names
        self.aux_comps = []
           
        # self.first_reading = True
        self.set_typical_geomag_range_filter() # ToDo: in case of IAGA format, set sepcial limits
                 
         
     def set_typical_geomag_range_filter(self):
        assert self.range_filter['min'] is None
        assert self.range_filter['max'] is None
        
        MIN_F_FIELD = 15000. # South Atlantic Anomaly core minimum is about 22 uT
        MAX_F_FIELD = 80000. # 88888.88 has a special measning in IAGA files, we should be below this
        
        if self.parent_mag_obs_data_set.instr_frame=="dIdD":
            MIN_VAL = MIN_F_FIELD
            MAX_VAL = MAX_F_FIELD * 1.2
        else:
            MIN_VAL = -MAX_F_FIELD
            MAX_VAL = MAX_F_FIELD
        
        self.range_filter['min'] = [ MIN_VAL ] * self.nof_comps
        self.range_filter['max'] = [ MAX_VAL ] * self.nof_comps
        if "F" in self.comp_names:
            F_index = self.comp_names.index("F")
            self.range_filter['min'][F_index]=MIN_F_FIELD
        
         
     def load_data(self, time_label_def, result_data_type = dat_file_reader.file_read_result_data_type.records_for_timestamps, time_shift_sec = None, in_case_of_missing_error = "none", return_array_format = "structured-array"):
        
        # Adding the new reader info to the list (and deleting the first element)
        #  _latest_used_dat_file_readers = _latest_used_dat_file_readers[1:] + [ [self.dat_file_def, _dat_file_reader]   ]

        # Reading in data file
        corrected_time_label_def = time_label_def if time_shift_sec is None else np.copy(time_label_def)+time_shift_sec/86400.0
        errmsg, timelabel_from_file, data_from_file, noncomplete_record_indeces, extra_data = self.get_data( corrected_time_label_def, result_data_type, in_case_of_missing_data = "none")
        
        # missing_record_indeces = []
        #if (isinstance(time_label_def, float)):
        #    if np.isnan(timelabel_from_file):
        #        missing_record_indeces = [0]
        #        missing_warn_msg = "{} is missing for time-label of {}".format(self.data_descr, time_label_def)
        #else:
        #    if (np.isnan(timelabel_from_file).any()):
        #        missing_record_indeces = np.nonzero( np.isnan(timelabel_from_file))[0]
        #        missing_warn_msg = "{:d} {} record{} missing for the total time-period  of {}({})...{}({}) ".format( len(missing_record_indeces), self.data_descr, ' is' if len(missing_record_indeces)==1 else 's are', mu.exceltime_2_datetimestr(time_label_def[0]), time_label_def[0], mu.exceltime_2_datetimestr(time_label_def[-1]), time_label_def[-1])

        if len(noncomplete_record_indeces)>0:
            missing_warn_msg = "{:d} {} record{} missing for the total time-period  of {}({})...{}({}) ".format( len(noncomplete_record_indeces), self.data_descr, mu.plural_suffix(len(noncomplete_record_indeces),'is','s are'), mu.exceltime_2_datetimestr(time_label_def[0]), time_label_def[0], mu.exceltime_2_datetimestr(time_label_def[-1]), time_label_def[-1])    
        
            if in_case_of_missing_error != "none":
                assert np.all( np.isnan( data_from_file[noncomplete_record_indeces[0]] )  )
                pass
            else:
                print(  "{} in directory \"{:}\"".format( missing_warn_msg, self.path ) )
                if return_array_format=="structured-array":
                    assert len(time_label_def) == len(timelabel_from_file) # The listing below is based on the assumption that 'timelabel_from_file' corresponds to "time_label_def"
                    mu.list_timelabels("Missing records:",  time_label_def[np.isnan(timelabel_from_file)], 10   )

                if in_case_of_missing_error == "drop-error": assert False, "Fatal error!"

        assert len(self.mag_comps) + len(self.aux_comps) == data_from_file.shape[1]

        all_comps_list = []
        for comp in self.mag_comps: 
            assert not comp in all_comps_list
            all_comps_list += [ comp ]
        # if self.convert_to_instr_frame is not None:
        #     for new_instr_frame_comp in  self.convert_to_instr_frame:
        #         if not new_instr_frame_comp in all_comps_list:
        #             all_comps_list += [ new_instr_frame_comp ]
        # for calc_comp in calc_comps:
        #     if not calc_comp in all_comps_list: all_comps_list += [ calc_comp ]
        for aux_comp_name in self.aux_comps:
            assert not aux_comp_name in all_comps_list 
            all_comps_list += [ aux_comp_name ]
               
        if return_array_format=="structured-array":
            # User requested more records - we create a structured array
            dtype = np.dtype( [(comp, float) for comp in all_comps_list] )
            result = np.empty( (len(timelabel_from_file)), dtype=dtype ) # to initialize with nan's: 'np.full( ... , fill_value = report_for_no_data)

            for comp_index, comp_name  in enumerate(self.mag_comps):
                result[comp_name][:] = data_from_file[:,comp_index]
            for aux_comp_name, comp_index in enumerate(self.aux_comps):
                result[aux_comp_name][:] = data_from_file[:,comp_index + len(self.mag_comps)]

        else:
            # Returning a single record as a dictionary
            result = dict()
            for comp_index, comp_name in enumerate(self.mag_comps):
                result[comp_name] = data_from_file[:,comp_index].copy()

        #if (not "F" in ref_existing_comps) and (not np.isnan(default_F_value)):
        #    ref_existing_comps += [ "F" ]
        #    mu.warning("{}data files does not contain F data, using default F value of {:.2f}".format(data_descr,default_F_value), True)
        #    if isinstance(ref_X,float):
        #        ref_F = default_F_value
        #    else:
        #        ref_F = np.empty([ len(ref_X) ])
        #        ref_F[:] = default_F_value

        # if self.baselines is not None:
        #     assert False
        #     assert (self.X0 is None) and (self.Y0 is None) and (self.Z0 is None) and (self.H0_nT is None) and (self.D0_nT is None)
        #     if self.frame == "XYZ":
        #         self.X0 = self.baselines[0]
        #         self.Y0 = self.baselines[1]
        #     elif self.frame == "HDZ":
        #         self.H0_nT = self.baselines[0]
        #         self.D0_nT = self.baselines[1]
        #     else:
        #         assert False
        #     self.Z0=self.baselines[2]

        #     self.baselines = None # We parsed/translated 'baselines', so making it invalied

        # Executing user's correction (is any)
        #if self.user_corr_func is not None:
        #    result = self.user_corr_func(result, self.file_mag_comps_str)
            
        # if self.convert_to_instr_frame is None:
        #     self.usage_instr_frame = self.file_instr_frame
        # else:
        #     if (self.convert_to_instr_frame == "XYZ") :
        #         if self.magnetometer_type == "DIDD":
        #             assert self.file_instr_frame == "dIdD"
        #             import sys
        #             sys.path.insert(1, 'p:/dev/datacq/_DIDD/')            
        #             import didd_utils
        #             didd_i0 = didd_utils.def_i0_get( time_label_def[0] )
        #             didd_d0 = didd_utils.def_d0_get( time_label_def[0] )
        #             sink_I, sink_D, didd_X, didd_Y, didd_Z = didd_utils.didd_process_deflected_values( time_label_def, result['F'], result['Dp'], result['Dm'], result['Ip'], result['Im'], didd_i0, didd_d0  )
        #             result['X']=didd_X
        #             result['Y']=didd_Y
        #             result['Z']=didd_Z
        #             self.instr_frame = "XYZ"
        #             # all_comps_list += ['X','Y','Z']
        #         elif self.file_instr_frame in ["HDZ", "HDZF"]:
        #             assert self.parent_mag_obs_data_set is not None
        #             assert self.parent_mag_obs_data_set.file_decbas is not None
        #             result['X'], result['Y'] = gu.HD_2_XY(result['H'], result['D'], self.parent_mag_obs_data_set.file_decbas)
        #             # After conversion, decbas will be zero:
        #             self.parent_mag_obs_data_set.usage_decbas = 0.0
        #         else:
        #             assert False, "Feature not supported"
        #     else:
        #         assert False, "Feature not supported"
        #     self.usage_instr_frame = self.convert_to_instr_frame
        #     # 2026.02: Instr. frame does not have 'F' if (not "F" in  self.usage_instr_frame) and ("F" in self.file_instr_frame): self.usage_instr_frame += "F"
            

        # Auto-calculating missing components
        # First, 'H' must be calcualted, because it is needed for some other calcualtions:

        # Applying known scale and baseline values before component conversion.
        # It is required for the calculation of Inclination/Declination angles

        # if self.A_scale is not None:
        #     result[ self.instr_frame[0] ] *=  self.A_scale
        #     # print( self.A_scale )
        # if self.B_scale is not None:  result[ self.instr_frame[1] ] *=  self.B_scale
        # if self.Z_scale is not None:  result[ 'Z' ] *=  self.Z_scale

        # assert not hasattr(self, 'H0'),    "Use 'baselines' dictionary"
        # assert not hasattr(self, 'H0_nT'), "Use 'baselines' dictionary (H0_nT={})".format(self.H0_nT)
        # assert not hasattr(self, 'D0'),    "Use 'baselines' dictionary"
        # assert not hasattr(self, 'D0_nT'), "Use 'baselines' dictionary"
        # assert not hasattr(self, 'Z0'), "Use 'baselines' dictionary"
        # assert not hasattr(self, 'X0'), "Use 'baselines' dictionary"
        # assert not hasattr(self, 'Y0'), "Use 'baselines' dictionary"
        # assert not hasattr(self, 'decbas'), "Use 'baselines[COMPNAME__BASELINE__D_ANGLE]' dictionary item or file_decbas!" 

        # self.applied_baselines_at_file_readings = {}
        # if self.baselines is not None:
        #     if gu.COMPNAME__BASELINE__HX_SENSOR in self.baselines.keys():
        #         result[ gu.COMPNAME_VARIOMETER_H ] +=  self.baselines[gu.COMPNAME__BASELINE__HX_SENSOR]
        #         self.applied_baselines_at_file_readings[gu.COMPNAME__BASELINE__HX_SENSOR] = self.baselines[gu.COMPNAME__BASELINE__HX_SENSOR]
        #     # if self.Y0 is not None:
        #     #     result[ 'Y' ] +=  self.Y0
        #     #     self.Y0_applied_at_file_reading = self.Y0
        #     # else:
        #     #     self.Y0_applied_at_file_reading = None
        #     if gu.COMPNAME__BASELINE__DY_SENSOR in self.baselines.keys():
        #         assert False, "Feature not implemented"

        #     if gu.COMPNAME__BASELINE__Z_SENSOR in self.baselines.keys():
        #         result[ gu.COMPNAME_VARIOMETER_Z ] +=  self.baselines[gu.COMPNAME__BASELINE__Z_SENSOR]
        #         self.applied_baselines_at_file_readings[gu.COMPNAME__BASELINE__Z_SENSOR] = self.baselines[gu.COMPNAME__BASELINE__Z_SENSOR]
        #     # if self.Z0 is not None:
        #     #     result[ 'Z' ] +=  self.Z0
        #     #     self.Z0_applied_at_file_reading = self.Z0
        #     # else:
        #     #     self.Z0_applied_at_file_reading = None

        #     # if self.H0_nT is not None:
        #     #     result[ gu.COMPNAME_VARIOMETER_H ] +=  self.H0_nT
        #     #     self.H0_nT_applied_at_file_reading = self.H0_nT
        #     # else:
        #     #     self.H0_nT_applied_at_file_reading = None

        #     # if self.D0_nT is not None:
        #     #     result[ gu.COMPNAME_VARIOMETER_D ] +=  self.D0_nT
        #     #     self.D0_nT_applied_at_file_reading = self.D0_nT
        #     # else:
        #     #     self.D0_nT_applied_at_file_reading = None



        # # if self.F_pillar_diff is not None:
        # #     assert 'F' in file_mag_comps
        # #     result['F'] += self.F_pillar_diff # F_abs = F_variometer + F_pillar_diff
        # #     self.F_pillar_diff_applied_at_file_reading = self.F_pillar_diff
        # # else:
        # #     self.F_pillar_diff_applied_at_file_reading = None

        # if gu.COMPNAME_TOTAL_H in calc_comps:
        #     if self.instr_frame[:4] == "XYZ":
        #         result[gu.COMPNAME_TOTAL_H] = gu.XY_2_H( result['X'], result['Y'] )
        #     elif self.instr_frame[:4] == "HDZ":
        #         result[gu.COMPNAME_TOTAL_H] = gu.XY_2_H( result[gu.COMPNAME_VARIOMETER_H], result[gu.COMPNAME_VARIOMETER_D] )
        #     else:
        #         result[gu.COMPNAME_TOTAL_H] = np.nan
        #         assert False

        # for calc_comp in calc_comps:
        #     if not calc_comp in self.mag_comps:
        #         if (calc_comp==gu.COMPNAME_TOTAL_H):
        #             continue # It was already calcualted, see above
        #         elif (calc_comp==gu.COMPNAME_DECL_ANGLE):
        #             if self.usage_instr_frame == "XYZ":
        #                 result[gu.COMPNAME_DECL_ANGLE] = gu.XY_2_Decl( result['X'], result['Y'] )
        #             elif self.usage_instr_frame == "HDZ":
        #                 result[gu.COMPNAME_DECL_ANGLE] = gu.HE_2_Decl( result[gu.COMPNAME_VARIOMETER_H], result[gu.COMPNAME_VARIOMETER_D],  self.baselines[gu.COMPNAME__BASELINE__D_ANGLE] )
        #             else:
        #                 assert False
        #         elif (calc_comp==gu.COMPNAME_INCL_ANGLE):
        #             if self.usage_instr_frame == "XYZ":
        #                 if 'F' in self.mag_comps:
        #                     result[gu.COMPNAME_INCL_ANGLE] = np.arcsin(result['Z'] / result['F'])
        #                 else:
        #                     result[gu.COMPNAME_INCL_ANGLE] = gu.totalH_Z__2__Incl( result[gu.COMPNAME_TOTAL_H], result['Z'] )
        #             elif self.usage_instr_frame == "HDZ":
        #                 result[gu.COMPNAME_INCL_ANGLE] = gu.totalH_Z__2__Incl( result[gu.COMPNAME_TOTAL_H], result['Z'] )
        #             else:
        #                 result[gu.COMPNAME_INCL_ANGLE] = np.nan
        #                 assert False
        #         elif (calc_comp=="X"):
        #             if self.usage_instr_frame == "HDZ":
        #                 assert False 
        #                 assert (self.usage_decbas is not None)
        #                 result['X'] = result[gu.COMPNAME_TOTAL_H] * np.cos(self.decbas)
        #             else:
        #                 assert False
        #         elif (calc_comp=="Y"):
        #             if self.usage_instr_frame == "HDZ":
        #                 assert False
        #                 assert (self.usage_decbas is not None)
        #                 result['Y'] = result[gu.COMPNAME_TOTAL_H] * np.sin(self.decbas)
        #             else:
        #                 assert False
        #         elif (calc_comp=="F"):
        #             result['F'] = np.sqrt(np.square(result[self.usage_instr_frame[0]]) + np.square(result[self.usage_instr_frame[1]]) + np.square(result[self.usage_instr_frame[2]]) )
        #         else:
        #             assert False, "Unknown calculated component: {}".format(calc_comp)

        # #if isinstance(self.mag_comps, str):
        # #    all_comps = self.file_mag_comps + calc_comps
        # #else:
        # #    for calc_comp in calc_comps:
        # #        all_comps += [calc_comp]
                

        # if False:
        #     if self.decbas is not None:
        #         assert False, "Use 'baselines[mu.COMPNAME__BASELINE__D_ANGLE]' instead"
        #         if 'X' in result.dtype.names:
        #             result[ 'X' ] =  result[ 'X' ] * np.cos(self.decbas)
        #             result[ 'Y' ] =  result[ 'Y' ] * np.cos(self.decbas)
        #         else:
        #             result[ 'H' ] =  result[ 'H' ] * np.cos(self.decbas)
        #             result[ 'D' ] =  result[ 'D' ] * np.cos(self.decbas)


        # if hasattr(self, "Incl0"):
        #     assert False
        #     if self.Incl0 is not None:
        #         assert np.isfinite(self.Incl0)
        #         if 'I' in result.dtype.names:
        #             result['I'] = result['I']  + self.Incl0
        #         else:
        #             mu.warning("Can not apply I pillar difference of {:.3f} degree for {}".format( math.degrees(self.Incl0), data_descr ), True)
        # if hasattr(self, "Decl0"):
        #     assert False
        #     if self.Decl0 is not None:
        #         assert np.isfinite(self.Decl0)
        #         if 'D' in result.dtype.names:
        #             result['D'] = result['D']  + self.Decl0
        #         else:
        #             mu.warning("Can not apply D pillar difference of {:.3f} degree for {}".format( math.degrees(self.Decl0), data_descr ), True)

        # if save_filename is not None:
        #     mu.savetxt(save_filename,result)
        return result, noncomplete_record_indeces

         

# "Variometer + F" data set, or dIdD data set, of F only data set.
# Used as input for DIM baseline calculations
class  mag_obs_dataset:

    def __init__(self, dataset_def_id): # , variometer_aux_comps=[], data_path = None, F_data_is_from_DIM_data_file = False, aux_file_readers = []):

        # dataset_def_id = dataset_def_id.split("/") # Typically it is one id, but two id's can be concatenated: first idis the variometer, second part is the scalar file reader
        # Further items are for aux files
        # assert len(source_def_parts) in [1,2]
        #if data_path is None:
        #    data_path = [None]*len(source_def_parts)
        # assert len(data_path) == len(source_def_parts)
        # self.variometer_data_def = dataset_def_id[0]
        
        assert mag_obs_dataset_defs_yaml_path is not None, "Obserbatrory data sets configuration file is not yet processed"
        assert dataset_def_id in mag_obs_dataset_def_list.keys(), "No \"{}\" dataset found in file \"{}\"".format(dataset_def_id, mag_obs_dataset_defs_yaml_path)
        dataset_def = mag_obs_dataset_def_list[dataset_def_id]
        
        self.dataset_name = dataset_def["dataset-name"]
        
        self.instr_frame = dataset_def["variometer-orientation"]
        assert self.instr_frame in ["HDZ","XYZ","dIdD"]
        self.usage_frame = dataset_def["usage-frame"]
        assert self.usage_frame in ["HDZ","XYZ"]
        
        self.file_readers = []
        
        self.mag_comps_from_files = []
        self.aux_comps = []
        # self.F_data_source = F_data_source_enum.none
        
        self.usage_decbas = 0.0
        
        self.adjustments_dict = dataset_def.get("adjustments", {})
        # 'self.adjustments_dict' should be a dictionary
        assert type(self.adjustments_dict) == dict
        self.adjustment_dict_not_used_keys = set( self.adjustments_dict.keys() ) 
        interpolated_data_result_halfwindow_size = dataset_def.get("interpolation-window[sec]", None)
        if interpolated_data_result_halfwindow_size is not None: interpolated_data_result_halfwindow_size = interpolated_data_result_halfwindow_size / 86400.
                
        for data_file_def in dataset_def["data-files"]:
            
            file_reader = mag_obs_data_file_reader(data_file_def, self)
            file_reader.set_processing_parameters(raw_data_result_timestamp_diff = "auto", interpolated_data_result_halfwindow_size = interpolated_data_result_halfwindow_size ) 
            
            for magcomp in file_reader.mag_comps:
                assert not magcomp in self.mag_comps_from_files, f'Magnetic component "{magcomp}" is read from more sources'
                self.mag_comps_from_files  += [ magcomp ]
                #if magcomp == "F": 
                #    assert self.F_data_source == F_data_source_enum.none
                #    self.F_data_source = F_data_source_enum.obs_data_file
                
            # In case of aux comps, we allow repeating component names
            self.aux_comps += file_reader.aux_comps
            self.file_readers += [file_reader]
        
        self.mag_comps_calculated = []
        if (self.instr_frame != self.usage_frame):
            for calc_comp in self.usage_frame:
                if not calc_comp in self.mag_comps_from_files:
                    self.mag_comps_calculated += [ calc_comp ] 
                    
        self.mag_comps = self.mag_comps_from_files + self.mag_comps_calculated
        
        if self.instr_frame == "dIdD":
            import didd_utils
            self.DIDD_i0 = angle_utils.str_to_radians( dataset_def["dIdD-I0"] )
            self.DIDD_d0 = angle_utils.str_to_radians( dataset_def["dIdD-D0"] )
            self.DP_compname = didd_utils.get_deflected_component_name("D+", self.mag_comps_from_files) 
            self.DM_compname = didd_utils.get_deflected_component_name("D-", self.mag_comps_from_files) 
            self.IP_compname = didd_utils.get_deflected_component_name("I+", self.mag_comps_from_files) 
            self.IM_compname = didd_utils.get_deflected_component_name("I-", self.mag_comps_from_files) 
            
        # self.main_interval_def =  self.file_readers[0].file_record_interval_type



    def has_variometer_data(self):
        return True # We assume that is has variometer data. ToDo: implement

    def get_decbas(self):  
        return self.usage_decbas
    #    if self._decbas is not None: return self._decbas
    #    if self.variometer_data.instr_frame[:4] == "XYZ": return 0.0
    #    assert False  # no decbas data


    def has_F_data(self):
        return "F" in self.mag_comps # isinstance(self.F_data, mag_obs_data_reader) or (self.F_data ==  gu.F_DATA_SOURCES.VARIOMETER_DATA)
    
    def calc_F(self, load_data):
        if self.usage_frame == "XYZ":
            return gu.XYZ_2_F(load_data['X'], load_data['Y'], load_data['Z'] ) 
        else:
            assert False

    # def calc_frame(self):
    #    return self.variometer_data.calc_frame

    # Possible return values: XYZ, XYZF, HDZ, HDZF, and maybe: F, dIdD
    #def instr_frame(self):
    #   result = self.variometer_data.usage_instr_frame
    #   #if isinstance(self.F_data, mag_obs_data_reader):   2026.02: 'instr_framew' does not have 'F'
    #   #    assert not "F" in result
    #   #    result = result + "F"
    #   return result
    
    
    def load_data(self, b_load_mag_data, b_load_aux_data, time_label_def, meas_pos = None, report_for_no_data = np.nan, in_case_of_missing_error = "none", save_filename = None, time_shift_sec = None, force_using_time_interpolation = False, info_lines_prefix = "", result_data_type = dat_file_reader.file_read_result_data_type.records_for_timestamps):
        
        result_dat = None
        result_noncomplete_record_indeces = None
        
        if isinstance(time_label_def, float):
            return_array_format = "dictionary"  # One record as a dictionaly
        else:
            return_array_format = "structured-array"
            

        if (len(self.aux_comps)==0) or  (b_load_mag_data and not b_load_aux_data): loaded_data_type =  "magnetic" 
        elif (not b_load_mag_data and b_load_aux_data): loaded_data_type =  "auxiliary" 
        else: loaded_data_type =  "all types of " 
        
        read_from_file_readers = []
        for file_reader in self.file_readers:
            if (b_load_mag_data and  (len(file_reader.mag_comps)>0)) or ( b_load_aux_data and  (len(file_reader.aux_comps)>0)  ):
                read_from_file_readers += [ file_reader]
                
                
        assert len(read_from_file_readers) > 0
        
        intro_line = "Loading {} data from \"{}\" dataset".format(loaded_data_type, self.dataset_name)
        
        self.info_lines_prefix = info_lines_prefix
        file_reader.info_lines_prefix = info_lines_prefix
        file_reader.detailed_file_open_line = False
        
        if info_lines_prefix == "":
            print("----- {} -----".format(intro_line))
            step_header = "-- "
        else:
            print(info_lines_prefix+intro_line)
            step_header = info_lines_prefix+"  --"
            
        for index, file_reader in enumerate(read_from_file_readers):
            read_comps = []
            if b_load_mag_data: read_comps += file_reader.mag_comps
            if b_load_aux_data: read_comps += file_reader.axu_comps
            assert len(read_comps) > 0
            if len(read_from_file_readers) > 1: 
                print(step_header + "{}/{} step: reading files to get {} components".format(index+1,len(read_from_file_readers), read_comps ) )
                
            if (time_shift_sec is None) and (not force_using_time_interpolation) and (file_reader.timeindex_interval_in_sec==1):
                # Reading 1 second data 
                file_reader.set_processing_parameters(raw_data_result_timestamp_diff = "auto", interpolated_data_result_halfwindow_size = None) 
                # file_read_method = dat_file_reader.file_read_method_enum.exact_time_match
            else:
                # Reading 1 minute files or 5 secpnd files
                # assert file_reader.timeindex_interval_in_sec==60
                file_reader.set_processing_parameters(raw_data_result_timestamp_diff =  "auto", interpolated_data_result_halfwindow_size = 61.0 / 86400.) 
                # file_read_method = dat_file_reader.file_read_method_enum.raw_interpolated
                
            result_dat, result_noncomplete_record_indeces = file_reader.load_data(time_label_def, result_data_type = result_data_type, time_shift_sec = time_shift_sec, in_case_of_missing_error = in_case_of_missing_error, return_array_format = return_array_format)
               
            assert len(read_from_file_readers) == 1, "Reading more different files feature is temporarly unavailable+"
         

        result_dat, result_noncomplete_record_indeces = self.calculate_derived_magnetic_components(time_label_def, result_dat, result_noncomplete_record_indeces)
        
        # Applying 'self.adjustments_dict' after loading:
        for key, value in self.adjustments_dict.items():
            if key[-7:] == "-offset": # " {Z-offset: -1.0} "
                comp = key[:-7]
                result_dat[comp] += value
                self.adjustment_dict_not_used_keys.discard(key)
                        
                        
        return result_dat, result_noncomplete_record_indeces
    
    
    def calculate_derived_magnetic_components(self, timestamps, data, noncomplete_record_indeces):
        
        if len(self.mag_comps_calculated) == 0: return data, noncomplete_record_indeces
        
        if self.instr_frame == "dIdD":
            import didd_utils
            if self.usage_frame == "XYZ":
                _, _, X,Y,Z  = didd_utils.didd_process_deflected_values( timestamps, data["F"], data[self.DP_compname], data[self.DM_compname], data[self.IP_compname], data[self.IM_compname], self.DIDD_i0, self.DIDD_d0 )
                data = rfn.append_fields(data, ['X','Y','Z'], data=[X,Y,Z], dtypes =  ("f8", "f8", "f8"), usemask = False )
                # Resoulting data may have new nan's results, where DIDD calculation failes!                
                new_noncomplete_record_indeces = np.where(np.isnan(X))
                new_noncomplete_record_indeces = np.union1d(noncomplete_record_indeces, new_noncomplete_record_indeces)
                
                
            elif self.usage_frame == "HDZ":
                pass
            else:
                assert False, f'Not supported frame: {self.calc_frame}'
            
        else:
            assert False, "Feature temporarly unavailable"
            
        return data, new_noncomplete_record_indeces


    def get_D0_nT_applied_at_file_reading(self):
        return 0.0 # self.variometer_data.get_D0_nT_applied_at_file_reading()

    def get_H0_nT_applied_at_file_reading(self):
        return 0.0 # self.variometer_data.get_H0_nT_applied_at_file_reading()

    def get_Z0_nT_applied_at_file_reading(self):
        return 0.0 # self.variometer_data.get_Z0_nT_applied_at_file_reading()


    def get_warnings(self):
        result = []
        if len(self.adjustment_dict_not_used_keys)>0:
            warn_msg = "Not used adjustment key{} in observatory data set: {}".format(mu.plural_suffix(len(self.adjustment_dict_not_used_keys)), mu.list_to_string(list(self.adjustment_dict_not_used_keys)))
            result += [ warn_msg ]
        return result
                
                 
    # Can be used to represent the variation in a time wondows
    # It tries to automatically find the optimal step
    def plot_data_for_date_range(self, date_from, date_to, mark_datetimes = None, datetime_separators = None):
        N_days = date_to - date_from
        if (N_days<1.0):
            timelabels = np.arange(date_from, date_to, 1.0/86400.0)
            result_data_type = dat_file_reader.file_read_result_data_type.records_for_timestamps
        else:
            timelabels = dat_file_reader.get_minute_timelabels(date_from, date_to)
            result_data_type = dat_file_reader.file_read_result_data_type._1min_arithm_avg
            
        mag_data_for_plot, missing_records = self.load_data( True, False, timelabels, result_data_type = result_data_type )
        self.plot_load_data(timelabels, mag_data_for_plot, "Data for entire period", mark_datetimes = mark_datetimes, datetime_separators = datetime_separators)
            

    # 'mark_datetimes'  can be a list of timelabels, whos position will be marked on the plot
    def plot_load_data(self, timelabels, load_data, title, mark_datetimes = None, datetime_separators = None):
        plotted_comps = []
        for comp in "XHDYZF":
            if comp in load_data.dtype.names: plotted_comps += comp
        assert len(plotted_comps)>0

        # fig = plt.figure(figsize=(preferences.FIGURE_SIZE[0], preferences.FIGURE_SIZE[1]*3), layout = ['constrained', 'compressed', 'tight', 'none', None][0] )
        num_plots = len(plotted_comps)
        fig, axes = plt.subplots(num_plots, 1, sharex=True, figsize=(preferences.FIGURE_SIZE[0], preferences.FIGURE_SIZE[1]*3))
        fig.suptitle(title)
        fig.subplots_adjust(hspace=0)

        for icomp, comp in enumerate(plotted_comps):
            # plt.subplot( len(plotted_comps), 1, icomp+1 )
            fig.sca( axes[icomp] )
            comp_color = preferences.magnetic_comp_plot_colors[comp]
            # b_configure_datetime_axis = (icomp==len(plotted_comps)-1) and (mark_datetimes is None)
            plot_utils.plot_datetime_y( timelabels,  [ load_data[comp] ], comp + " (nT)"  ,['*','-'][1], color=comp_color, b_close_plot=False )
            if mark_datetimes is not None:
                Y_closest = np.empty(len(mark_datetimes))
                for i,marked_datetime in enumerate(mark_datetimes):
                    assert marked_datetime>=timelabels[0],  "{} < {}".format(marked_datetime,timelabels[0])
                    assert marked_datetime<=timelabels[-1]+1/86400.0, "{} > {} ({})".format(marked_datetime,timelabels[-1], marked_datetime - timelabels[-1] )
                    closest_index =  (np.abs(timelabels - marked_datetime)).argmin()
                    Y_closest[i] = load_data[comp][closest_index]
                if comp_color == "red":
                    comp_color = "green"
                else:
                    comp_color = "red"
                plot_utils.plot_datetime_y( mark_datetimes,  [ Y_closest ], comp + " (nT)"  ,'X', color=comp_color, b_close_plot=False, vertical_separators = datetime_separators )
        plt.show()



def read_obs_datasets(i_mag_obs_dataset_defs_yaml_path):
    
    global mag_obs_dataset_defs_yaml_path
    mag_obs_dataset_defs_yaml_path = i_mag_obs_dataset_defs_yaml_path
    with open(mag_obs_dataset_defs_yaml_path, "r", encoding="utf-8") as _f:
        mag_obs_dataset_defs_yaml = yaml.safe_load(_f)

    for _mag_obs_data_set_id, _data_set_def in mag_obs_dataset_defs_yaml.items():    
        _data_set_def = dict(_data_set_def)
        mag_obs_dataset_def_list[_mag_obs_data_set_id] = _data_set_def

