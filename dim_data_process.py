# -*- coding: utf-8 -*-
"""

DIM_data_process.py

Module to process DIM (including ADS) measurements.
It is typically used to calculate spot baseline values.
This is the heart of amdp software package.

@author: László Merényi 
"""


import math
import numpy as np
import sys
import tabulate
import matplotlib.pyplot as plt
import types
# NOTE: 'lmfit' is imported lazily, only inside the functions that actually
# use it. These are only reached for the "LM-*"/"GN"/"LM[GN]" calc methods.
import misc_utils as mu
import cli_utils as cu
import model_utils
import dim_utils  as du
import geomag_utils as gu
import dim_measurement_loader
import amdp_base_formulas as amdp
import trigonometric_utils as tu
import dat_file_reader
import dim_azimuth
import plot_utils
import dim_instr_model as dim_instr
import angle_utils
import geomag_obs_dataset
import copy
import geomag_utils_preferences as preferences
from enum import Enum
import re

calc_params_plot_symbols = [ '.',  'x', '<', '>', 'v', 'o' ]

F_DIFF_COLOR = "brown"
F_DIFF_FMT  = "x"


PRE_BASELINE_CALCULATION_VARIOMETER_CALIBRATION_SOURCES = ["none",  "user-defined", "calcvars-prev-results", "prev-calcvars-results" ]
PRE_BASELINE_CALCULATION_DIM_SENSOR_PARAM_SOURCES       = ["none", "calcvars-prev-results", "prev-calcvar-results" ]

CYCLE_FILTER__PILLAR_NAME = "pillar_name"

DEF_NR_OF_RANDOM_STEPS = 20

CALC_ERROR_CODE__NONE  = 0
CALC_ERROR_CODE__MISC_ERROR = 99
CALC_ERROR_CODE__EXCLUDED_BY_FILTER = 1
CALC_ERROR_CODE__UNEQUAL_POSITIONS  = 2
CALC_ERROR_CODE__SOLVERS_ERROR      = 3
CALC_ERROR_CODE__PILLAR_MISMATCH    = 4
CALC_ERROR_CODE__NO_VALID_REFMARK   = 5
CALC_ERROR_CODE__MISSING_VARDATA    = 6
CALC_ERROR_CODE__MISSING_FDATA      = 7
CALC_ERROR_CODE__TOO_BIG_RESID_STD  = 8
CALC_ERROR_CODE__TOO_FEW_M_DATA     = 9


CALC_ERROR_CODE__POST_CHECK_ERROR_CODES_START = 10
CALC_ERROR_CODE__POST_CHECK_1       = 10
CALC_ERROR_CODE__POST_CHECK_2       = 11
CALC_ERROR_CODE__POST_CHECK_I_1     = 12
CALC_ERROR_CODE__POST_CHECK_I_2     = 13
CALC_ERROR_CODE__POST_CHECK_H0      = 14

CALC_ERROR_CODE__NOT_CALCULATED     = 98

class F_data_source_enum(Enum):
     default = 0
     main_obs_dataset = 1
     separate_obs_dataset = 2 
     sites_default_F = 3
     DIM_data_file = 4
     calc_from_variometer = 5



_calc_params_last_id = -1

def create_copied_calc_var(base_calc_var):
    res = copy.copy(base_calc_var)
    global _calc_params_last_id
    _calc_params_last_id += 1
    res.id = str(_calc_params_last_id) # It must be a string, as it will be used as a differeniating component name
    res.result_file_name = base_calc_var.result_file_name +  "_{:02d}".format(_calc_params_last_id)
    return res



class spot_baselines_calculation_parameters_class:

    def __init__(self, params_dict, varcorr_method = "auto", meas_pos = None, DIM_instr = None, refmark_measured_data_sources = None, optimize_hc_eccentr_err = False,  optimize_vc_eccentr_err = False, mag_sensor_index = 0, nof_random_steps = DEF_NR_OF_RANDOM_STEPS, optimize_epsilon = True, epsilon_arcmin_startval = None, optimize_delta = True, delta_arcsec_startval = None, optimize_M_scale = True, M_scale_startval = None,  optimize_M_offset = True, M_offset_startval = None, mag_precorr_funcs = [], optimize_magcorr_funcs = [], cycle_filters= [], cycle_point_filters = [], test_mode = False, name = None, optimize_decbas = False, starting_mag_sensor = None, separate_D_I_processing = None): # , refpos_corr_data =  None
        # baseline_calc_approach = gu.BASELINE_CALCULATION_APPROACHES[0], 

        global _calc_params_last_id
        _calc_params_last_id += 1
        
        self.result_file_name = params_dict.pop("results-file", "amdp-results-{:02d}.dat".format(_calc_params_last_id))
        self.measpoints_result_file_name = params_dict.pop("pointwise-results-file", "amdp-measpoint-results-{:02d}.dat".format(_calc_params_last_id))        
        self.do_self_checks              = params_dict.pop("self-checks", True)
        F_source_str                     = params_dict.pop("F-source", "default") 
        
        
        calc_method = params_dict.pop("method")
        
        #self.baseline_calc_approach = baseline_calc_approach
        #assert self.baseline_calc_approach in gu.BASELINE_CALCULATION_APPROACHES
        self.decbas = None
        assert isinstance(meas_pos, gu.measurement_pos_class )
        
        if varcorr_method == "auto":
            if calc_method=="RESID":
                varcorr_method = "none"   # Lauridsen Residual method does not need variation correction
            elif calc_method in ["GN", "LM-DI"]:
                varcorr_method = ["avg","T[0]"][0]
            elif calc_method  == "LM/ID-var":
                varcorr_method = "LM/ID-var"
                assert False, "Check this"
            elif calc_method[:6] == "LM-BL/":
                varcorr_method = "none"
            else:
                assert False, f"Unknown calc-method: {varcorr_method}"

        if (name is None) or (name == ""):
            # name = calc_method
            name = "#{}".format(_calc_params_last_id)
            if varcorr_method!="none": name += "/w-varcorr"
        self.name = name
        self.id = str(_calc_params_last_id) # It must be a string, as it will be used as a differeniating component name
        self.method                = calc_method

        self.mag_sensor_index      = mag_sensor_index
        self.random_steps          = nof_random_steps
        self.optimize_M_scale      = optimize_M_scale
        self.M_scale_startval      = M_scale_startval
        self.optimize_M_offset     = optimize_M_offset
        self.M_offset_startval     = M_offset_startval
        self.optimize_epsilon      = optimize_epsilon
        self.optimize_delta        = optimize_delta
        self.epsilon_arcmin_startval= epsilon_arcmin_startval
        self.delta_arcsec_startval  = delta_arcsec_startval
        self.optimize_ortho_error   = False   # DI model can not optimize ortho error
        self.hc_precorr_funcs       = None
        self.vc_precorr_funcs       = None
        self.mag_precorr_funcs      = mag_precorr_funcs   # Fix (not optimized) mag sensor correction
        self.optimize_decbas        = optimize_decbas
        for magcorr_func in mag_precorr_funcs:
            if not isinstance(magcorr_func,dict):
                assert  callable( getattr(magcorr_func, "calc", None) )
        self.optimize_magcorr_funcs   = optimize_magcorr_funcs   # It can be empty, or one or more fixed correction and/or optimizable correction
        self.optimize_hc_eccentr_err  = optimize_hc_eccentr_err
        self.optimize_vc_eccentr_err  = optimize_vc_eccentr_err
        self.cycle_filters            = cycle_filters # Filters in/out cycles, e.g. with a pillar name
        self.cycle_point_filters      = cycle_point_filters # Filters out points from cycles
        self.opt_range = None  # Introduced in 2026.07.14, for LM method ranges. Format: dictionary, e.g. {"LM-A0:(-100,+100),  LM-B0:(-100,+100), LM-Z0:(40000,41000)}

        self.allow_unequal_measurement = False

        self.max_incl_change_deg = 3
        self.max_decl_change_deg = 3
        
        # self.min_nr_of_M_samples = 99999

        # self.refpos_measured_data_source__external_table_filename = _def_refpos_measured_data_source__external_table_filename

        if (calc_method in [ "NM", "RESID" ]) and (separate_D_I_processing is None):
            separate_D_I_processing = True # It gives better residuals in case of hc eccentr error
        self.separate_D_I_processing = separate_D_I_processing

        # assert gu.default_meas_pos is not None # The default measurement position must be selected prior calling this
        
        self.calc_residuals_from_last_result = False 
        self.calculated_null_positions_file = None
        aux_calculations = params_dict.pop("aux-calculations", [])
        for aux_calc_item in aux_calculations:
            if type(aux_calc_item) is dict:
                (aux_calc_item, aux_calc_param), = aux_calc_item.items()
                if aux_calc_item == "null-positions": 
                    self.calculated_null_positions_file = aux_calc_param
                else:
                    assert False
            else:
                if aux_calc_item=="residuals-from-last-window-results":
                    self.calc_residuals_from_last_result = True
                else:
                    assert False
        
        
        self.last_results  = { 'I':meas_pos.I_typical, 'D': meas_pos.D_typical, "baseline" : None }
        if starting_mag_sensor is not None:
            assert isinstance( starting_mag_sensor, dim_instr.mag_sensor_model_class)
            self.last_results["mag_sensor"] = copy.deepcopy(starting_mag_sensor)

        # Checking if varcorr_method is allowed for the calc_method
        if varcorr_method != "none":
            if calc_method in [ "LM-DI", "LM-D", "LM-|" ]:
                assert varcorr_method in  [ "LM/ID-var", "avg", "T[0]", "LM/dyn" ]
            elif calc_method in [ "NM", "RESID", "GN" ]:
                assert varcorr_method in  [ "none", "T[0]", "avg" ]
            elif calc_method in [  "LM[GN]" ]:
                assert varcorr_method in  [  "T[0]", "avg" ]  # At least one simple varcorr must be used
            else:
                assert False

        if self.separate_D_I_processing:
            assert calc_method in [ "NM", "RESID" ]

        self.varcorr_method = varcorr_method
        self.test_mode = test_mode

        self.refmark_measured_data_sources = refmark_measured_data_sources  # Defines the source(s) of reference mark reading data 
        
        self.true_angle_calculator = dim_azimuth.azimuth_calculator_def_class(meas_pos, dim_azimuth.pillar_data_source_enum.internal_pillar_data_table, references_sensor_pos = DIM_instr.hc_angle_reference_sensor_pos()  )

        if self.optimize_decbas: assert self.is_all_in_one_method()

        self.pre_baseline_calculation_variometer_calibration_source = PRE_BASELINE_CALCULATION_VARIOMETER_CALIBRATION_SOURCES[2]  # This can make the variation correction more accurate, in case of big epsilon/delta/scale-factor errors
        self.pre_baseline_calculation_dim_mag_sensor_parameters_source = PRE_BASELINE_CALCULATION_DIM_SENSOR_PARAM_SOURCES[0]

        
        if F_source_str[:15] == "mag-obs-dataset":   # Format of the setting: "mag-obs-dataset:MNY-DIDD-XYZF"
            F_source_str = F_source_str.replace(':',' ')
            dataset_id = F_source_str.split()[1]
            self.F_source_dataset = geomag_obs_dataset.mag_obs_dataset(dataset_id)
            F_source_str = "separate-mag-obs-dataset"            
        

        self.F_source =  {  \
            "default":                  F_data_source_enum.default,
            "main_obs_dataset" :        F_data_source_enum.main_obs_dataset,
            "separate-mag-obs-dataset": F_data_source_enum.separate_obs_dataset ,
            "sites-default-F":          F_data_source_enum.sites_default_F,
            "calc-from-variometer":     F_data_source_enum.calc_from_variometer
         } [F_source_str]
        
        self.baseline_tempco_corr = None
        self.baseline_tempco_corr_source = None

        # self.ae1_hc_cos_function = copy.deepcopy(ae1_hc_cos_function)
        # self.ae1_vc_cos_function = copy.deepcopy(ae1_vc_cos_function)
        # self.ae2_cos_function    = copy.deepcopy(ae2_cos_function)
        # if ae3_cos_function != None:
        #     assert ae3_cos_function.is_active()
        #     self.ae3_cos_function    = copy.deepcopy(ae3_cos_function)
        #     self.ae2_cos_function.omega_min = 1.2
        #     self.ae2_cos_function.omega_max = 1.9
        #     self.ae3_cos_function.omega_min = 1.91
        # else:
            #     self.ae3_cos_function  = None
            
        self.use_only_sample_indeces = params_dict.pop("use-only-samples", None)
        if self.use_only_sample_indeces is not None: 
            assert type(self.use_only_sample_indeces) == list
            
        precorr_funcs_def = params_dict.pop("angle-reader-calibration", None)
        if precorr_funcs_def:
            # "angle-reanding-corrections: {units: arcsec/deg, hc: 60/180 }"
            precorr_funcs_elements = [x.strip() for x in precorr_funcs_def.strip("'\" ").split(',')]
            
            for precorr_funcs_element in precorr_funcs_elements:
                if precorr_funcs_element.startswith("units:"):
                    assert precorr_funcs_element[6:] == "arcsec/deg"  # Unit of amplitude and pahsd
                    continue
                elif precorr_funcs_element.startswith("hc:"):
                    axis = "hc"
                elif precorr_funcs_element.startswith("vc:"):
                    axis = "vc"
                else:
                    assert False, "Invalid definition for angle calibration: " + precorr_funcs_element

                omega, amplitude_arcsec, phase_ = map(float, precorr_funcs_element[3:].split("/") )
                func =  tu.cos_function(  ampl = angle_utils.arcsec_to_rad(amplitude_arcsec), phase_ =  phase_, omega = omega )  
                if axis == "hc": 
                    if self.hc_precorr_funcs is None: self.hc_precorr_funcs = []
                    self.hc_precorr_funcs += [func]
                else:
                    if self.vc_precorr_funcs is None: self.vc_precorr_funcs = []
                    self.vc_precorr_funcs += [func]
                        
        self.angle_precorr_funcs = { "hc": self.hc_precorr_funcs, "vc" : self.vc_precorr_funcs }                                            
        
 
        self.opt_range = { }
        self.opt_start = { }
        
        opt_params = params_dict.get("parameters", [ ])
        assert type(opt_params) is list
        for opt_item in opt_params:
            assert type(opt_item) is dict
            assert len(list(opt_item.items() ))==1
            param_name = next(iter(opt_item))   
            param_def = opt_item[param_name]
            assert type(param_def) is dict
                
            fixed_value = param_def.pop("value", None)
            start_value = param_def.pop("start", None) #  It is optional
            assert not ( (start_value is not None) and (start_value is not None) )
            min_value = param_def.pop("min", None)
            max_value = param_def.pop("max", None)
            assert (min_value is None)  == (max_value is None)
            if (start_value is not None): assert min_value is not None
            if (start_value is not None): assert max_value is not None
            if (min_value is not None) and (start_value is None):
                start_value = (min_value + max_value) / 2.0
            
            match param_name:
                case "M-scale":
                    if fixed_value:
                        self.optimize_M_scale = False
                        assert np.isclose(fixed_value, 1.0)
                    else:
                        assert False
                case "H0" | "X0":                        
                    self.opt_start["BL-A0"] = start_value
                    self.opt_range["BL-A0"] = (min_value,max_value)
                case "D0" | "Y0":                        
                    self.opt_start["BL-B0"] = start_value
                    self.opt_range["BL-B0"] = (min_value,max_value)
                case "Z0":
                    self.opt_start["BL-Z0"] = start_value
                    self.opt_range["BL-Z0"] = (min_value,max_value)
                case _:
                    assert False
                        

        # Currently, calc_params_common has DIM time shift...
        # time_shift = params_dict.pop("DIM-time-shift"):
        #if time_shift:
        #    self.DIM_time_shift_sec
       

    def prepare_for_calc(self, mag_obs_dataset):
        
        if self.varcorr_method != "none":
            assert (mag_obs_dataset is not None)
            assert mag_obs_dataset.has_variometer_data() and mag_obs_dataset.has_F_data()
        #if (F_data_def  == F_DATA_IS_FROM_VARIOMETER_DATA ):
        #    assert (variometer_data_def is not None), "Variometer data is required for F values"

        #if self.baseline_calc_approach == "auto":
        #    if mag_obs_dataset is None:
        #        self.baseline_calc_approach = "none"
        #    elif mag_obs_dataset.instr_frame() in ["XYZ", "XYZF"]:
        #        self.baseline_calc_approach =  "Y0-nT"
        #    elif mag_obs_dataset.instr_frame() in ["HDZ", "HDZF"]:
        #        self.baseline_calc_approach =   "D0-angle"
        #    elif mag_obs_dataset.instr_frame() in ["dIdD"]:
        #        self.baseline_calc_approach =   "dIdD"
        #    else:
        #        assert False, mag_obs_dataset.instr_frame()
        #assert self.baseline_calc_approach in gu.BASELINE_CALCULATION_APPROACHES[1:], self.baseline_calc_approach


    def is_all_in_one_method(self):
        return (self.method[:6]=="LM-BL/") or ( self.method == "LM[GN]" )

    def name_str_get(self):
        return self.name
        # return "{:}*{:d}".format(self.method,  self.random_steps)


    def add_default_instr_model_params(self, params, instr):  # ~set_model_params()~
        import lmfit

        assert isinstance(params, lmfit.Parameters)
        mag_sensor = instr.mag_sensors[self.mag_sensor_index]
        M_scale = self.M_scale_startval if self.M_scale_startval is not None else mag_sensor.scale
        assert (M_scale > mag_sensor.scale_min) and (M_scale < mag_sensor.scale_max), f"{M_scale=} M_scale_min:{mag_sensor.scale_min} M_scale_max:{mag_sensor.scale_max} using-startval:{self.M_scale_startval}"
        # model.set_param_hint('M_scale',       value = M_scale,      min = mag_sensor.scale_min,     max = mag_sensor.scale_max,  vary = self.optimize_M_scale )     # Note: M_scale range should not cross zero, as everything may be inverted after
        params.add_many( lmfit.Parameter('M_scale',       value = M_scale,      min = mag_sensor.scale_min,     max = mag_sensor.scale_max,  vary = self.optimize_M_scale, user_data= ("nT/nT",)) )

        M_offset = self.M_offset_startval if self.M_offset_startval is not None else mag_sensor.offset
        assert (M_offset > mag_sensor.offset_min) and (M_offset < mag_sensor.offset_max)
        # model.set_param_hint('M_offset',     value = M_offset,     min = mag_sensor.offset_min,    max = mag_sensor.offset_max, vary = self.optimize_M_offset)
        params.add_many( lmfit.Parameter(        'M_offset',     value = M_offset,     min = mag_sensor.offset_min,    max = mag_sensor.offset_max, vary = self.optimize_M_offset, user_data= ("nT",)) )

        epsilon = angle_utils.arcmin_to_rad(self.epsilon_arcmin_startval) if self.epsilon_arcmin_startval is not None else mag_sensor.epsilon()
        assert (epsilon > angle_utils.arcmin_to_rad(mag_sensor.epsilon_min_arcmin)) and (epsilon < angle_utils.arcmin_to_rad(mag_sensor.epsilon_max_arcmin))
        # model.set_param_hint('epsilon',     value = epsilon,    min = angle_utils.arcmin_to_rad(mag_sensor.epsilon_min_arcmin),max = angle_utils.arcmin_to_rad(mag_sensor.epsilon_max_arcmin),   vary = self.optimize_epsilon )
        params.add_many( lmfit.Parameter(        'epsilon',     value = epsilon,    min = angle_utils.arcmin_to_rad(mag_sensor.epsilon_min_arcmin),max = angle_utils.arcmin_to_rad(mag_sensor.epsilon_max_arcmin),   vary = self.optimize_epsilon, user_data= ("'",  angle_utils.RAD_2_ARCMIN ) ) )

        delta = angle_utils.arcsec_to_rad(self.delta_arcsec_startval) if self.delta_arcsec_startval is not None else mag_sensor.delta()
        assert (delta > angle_utils.arcsec_to_rad(mag_sensor.delta_min_arcsec)) and (delta < angle_utils.arcsec_to_rad(mag_sensor.delta_max_arcsec))
        # model.set_param_hint('delta',       value = delta,      min = angle_utils.arcsec_to_rad(mag_sensor.delta_min_arcsec),  max = angle_utils.arcsec_to_rad(mag_sensor.delta_max_arcsec),     vary = self.optimize_delta)
        params.add_many( lmfit.Parameter('delta',       value = delta,      min = angle_utils.arcsec_to_rad(mag_sensor.delta_min_arcsec),  max = angle_utils.arcsec_to_rad(mag_sensor.delta_max_arcsec),     vary = self.optimize_delta, user_data=("\"", angle_utils.RAD_2_ARCSEC)) )

        self.add_eccentr_model_params(params, instr)

        # model.set_param_hint('ortho_error', value = angle_utils.arcsec_to_rad(instr.ortho_error_arcsec), min = -angle_utils.arcsec_to_rad(instr.ortho_error_max_arcsec),    max = angle_utils.arcsec_to_rad(instr.ortho_error_max_arcsec),  vary = self.optimize_ortho_error)
        params.add_many( lmfit.Parameter('ortho_error', value = angle_utils.arcsec_to_rad(instr.ortho_error_arcsec), min = -angle_utils.arcsec_to_rad(instr.ortho_error_max_arcsec),    max = angle_utils.arcsec_to_rad(instr.ortho_error_max_arcsec),  vary = self.optimize_ortho_error) )


    def add_eccentr_model_params(self, params, instr):
        import lmfit
        eccentr_ampl_min = 0.0
        hc_eccentr_ampl_startval = instr.hc_axis.eccentr_error_ampl
        vc_eccentr_ampl_startval = instr.vc_axis.eccentr_error_ampl
        if (self.optimize_hc_eccentr_err) and (np.isnan(instr.hc.eccentr_error_ampl)): hc_eccentr_ampl_startval = angle_utils.arcsec_to_rad(30)
        if (self.optimize_vc_eccentr_err) and (np.isnan(instr.vc.eccentr_error_ampl)): vc_eccentr_ampl_startval = angle_utils.arcsec_to_rad(30)

        # model.set_param_hint('hc_ecc_ampl',  value = hc_eccentr_ampl_startval,  min = eccentr_ampl_min, max = instr.hc.max_eccentr_error_ampl, vary = self.optimize_hc_eccentr_err)
        # model.set_param_hint('vc_ecc_ampl',  value = vc_eccentr_ampl_startval,  min = eccentr_ampl_min, max = instr.vc.max_eccentr_error_ampl, vary = self.optimize_vc_eccentr_err)
        params.add_many( lmfit.Parameter('hc_ecc_ampl',  value = hc_eccentr_ampl_startval,  min = eccentr_ampl_min, max = instr.hc_axis.max_eccentr_error_ampl, vary = self.optimize_hc_eccentr_err, user_data=model_utils.conv_radian_to_arcsec) )
        params.add_many( lmfit.Parameter('vc_ecc_ampl',  value = vc_eccentr_ampl_startval,  min = eccentr_ampl_min, max = instr.vc_axis.max_eccentr_error_ampl, vary = self.optimize_vc_eccentr_err, user_data=model_utils.conv_radian_to_arcsec) )

        # model.set_param_hint('hc_ecc_phase', value = instr.hc.eccentr_error_phase, min = 0.0, max = math.pi * 2.0 - 1e-6, vary = self.optimize_hc_eccentr_err)
        # model.set_param_hint('vc_ecc_phase', value = instr.vc.eccentr_error_phase, min = 0.0, max = math.pi * 2.0 - 1e-6, vary = self.optimize_vc_eccentr_err)
        params.add_many( lmfit.Parameter('hc_ecc_phase', value = instr.hc_axis.eccentr_error_phase, min = 0.0, max = math.pi * 2.0 - 1e-6, vary = self.optimize_hc_eccentr_err, user_data=model_utils.conv_radian_to_deg ) )
        params.add_many( lmfit.Parameter('vc_ecc_phase', value = instr.vc_axis.eccentr_error_phase, min = 0.0, max = math.pi * 2.0 - 1e-6, vary = self.optimize_vc_eccentr_err, user_data=model_utils.conv_radian_to_deg ) )

        # ECCENTR_ERR_AMPL__MAX   = angle_utils.arcsec_to_rad(  900  )
        # ECCENTR_ERR_AMPL__START = angle_utils.arcsec_to_rad(  60  )


    #
    # Can be called after an LM optimization
    def apply_instr_model_params(self, optimzied_params, instr):
        mag_sensor = instr.mag_sensors[self.mag_sensor_index]
        mag_sensor.scale = optimzied_params['M_scale'].value
        mag_sensor.offset = optimzied_params['M_offset'].value
        mag_sensor.epsilon_arcmin = math.degrees(optimzied_params['epsilon'].value) * 60
        mag_sensor.delta_arcsec = math.degrees(optimzied_params['delta'].value) * 3600.0

# Absolute measurements processing results are stored in 'result_array' list: each calculation variation will have a 'datarecord_array_class' storage
# This function calculates some information required to create and use this structured array.
def get_amdp_results_array_data(baseline_calc_approach, XH_compname, YD_compname, used_pillars, temp_comps, calc_params = None ):
    
    result_array_dtypes = [  ]    
    concatenated_plot_y_components = [ "resid_std" ]        
    result_array_comps_params = {  } 
    max_datetime_str_len = len(mu.exceltime_2_datetimestr(49000.0))
    timelabel_decimal_digits = mu.TIMELABEL_DECIMAL_DIGIT_FOR_SECOND # Deprecated: preferences.TIMELABEL_PREC    
    max_pillar_name_len = max( len(s) for s in used_pillars)
    
    def add_result_comp(add_to_table, comp_name, comp_type, comp_unit = "", decimal_digits = 0, optimized_param = False, plot = False):
        
        nonlocal result_array_dtypes, concatenated_plot_y_components, result_array_comps_params
        if decimal_digits is None: 
            if (comp_type[0]=='f') and (comp_unit!=""):
                decimal_digits = gu.get_decimal_digit_for_unit(comp_unit)
            else:
                decimal_digits = 0
        
        result_array_dtypes += [ (comp_name,  comp_type) ]   # To add comp_unit as metadata, use: np.dtype(comp_type, metadata={'unit': comp_unit} )
        # Old method: ra_cycles.comp_add(comp_name, comp_unit, decimal_digits, comp_type)        
        
        result_array_comps_params[comp_name] = types.SimpleNamespace()
        result_array_comps_params[comp_name].decimal_digits = decimal_digits
        result_array_comps_params[comp_name].add_to_table = add_to_table        
        result_array_comps_params[comp_name].unit = comp_unit
        result_array_comps_params[comp_name].optimized_param = optimized_param
        if comp_type[0]=='u':
            result_array_comps_params[comp_name].empty_value = 0xFFFFFFFF
            result_array_comps_params[comp_name].format_mask  = "%u"
        elif comp_type=='f4':
            result_array_comps_params[comp_name].empty_value = np.float32('nan')
            result_array_comps_params[comp_name].format_mask  = f"%.{decimal_digits}f"
        elif comp_type=='f8':
            result_array_comps_params[comp_name].empty_value = np.float64('nan')
            result_array_comps_params[comp_name].format_mask  = f"%.{decimal_digits}f"
        elif comp_type[0]=='U':
            result_array_comps_params[comp_name].empty_value = "na"
            result_array_comps_params[comp_name].format_mask  = "%s"
        else:
            assert False, f"{comp_name=} {comp_type=}"
        
        if isinstance(plot,bool):
            # 'plot' refers to a simple plot
            if plot and (not comp_name in concatenated_plot_y_components):
                concatenated_plot_y_components += [comp_name]
        else:
            # 'plot' should be a list of component names, that are plotted together like ["I","D"]
            assert isinstance(plot,list)
            concatenated_plot_y_components += [plot]
            
    add_result_comp(False, "win_nr",             "u4")   
    add_result_comp(False, "error_code",         "u4")
    add_result_comp(False, "program_type_idx",   "u4") # result_array[calc_params_index].record_data_set("main_program_type_idx", processing_window_index, meas_method_info[0])
    # add_result_comp(False, True, "win_name",           f"U{max_windown_name_len}"  )
    add_result_comp(False, "pier",               f"U{max_pillar_name_len}" )
    add_result_comp(False, "input_file_win_idx",  "u4" ) #  result_array[calc_params_index].record_data_set("file_cycle_id",  processing_window_index,  cycle_file_id_info[0] )
    add_result_comp(False, "datetime-start-str", f"U{max_datetime_str_len}")
    add_result_comp(False, "datetime-start",     "f8", decimal_digits = timelabel_decimal_digits)
    add_result_comp(False, "datetime-end",       "f8", decimal_digits = timelabel_decimal_digits)
    add_result_comp(False, "datetime-mid",       "f8", decimal_digits = timelabel_decimal_digits)
    add_result_comp(False, "datetime-Incl",      "f8", decimal_digits = timelabel_decimal_digits)
    add_result_comp(False, "datetime-Decl",      "f8", decimal_digits = timelabel_decimal_digits)
    # add_result_comp(False, True, "hours-of-day",       "f4", "h",  decimal_digits = 3)
    
    add_result_comp(True, "residual-stddev",            "f4", "nT",  decimal_digits = 1, plot = True) 
    add_result_comp(True, "azimark-hc_window-avg","f8", angle_utils.user_angle_unit, decimal_digits = angle_utils.user_angle_unit_decimal_digits  )
    add_result_comp(True, "azimark-hc-change_mean-removed",    "f8", "arcsec", decimal_digits = 1)
  
    
    add_result_comp(False, "duration",     "f4",  "min", decimal_digits = 1, plot = True ) 
    
    
    if baseline_calc_approach == "D0-angle":
        add_result_comp(True, "H0",                "f8", "nT",     decimal_digits = 2, plot = True )
        add_result_comp(True, "D0_nT_from_angle",  "f8", "nT",     decimal_digits = 2, plot = True )
        add_result_comp(True, "D0-angle",          "f8", "arcmin", decimal_digits = 2, plot = True  )
        baseline_components = [ "H0", "D0_nT_from_angle", "D0-angle" ]
    elif baseline_calc_approach == "Y0-nT":
        add_result_comp(True, XH_compname + "0",   "f8", "nT",  decimal_digits = 2, plot = True )
        add_result_comp(True, YD_compname + "0",   "f8", "nT",  decimal_digits = 2, plot = True )        
        add_result_comp(True, "D0-angle",          "f8", "min", decimal_digits = 2 , plot = True )        
        baseline_components = [ XH_compname + "0", YD_compname + "0", "D0-angle" ]
    elif baseline_calc_approach == "dIdD":  
        assert mag_obs_dataset.instr_frame() == "dIdD"
        assert False, "Not yet implemented"
    else:
        assert False
        
    add_result_comp(True, "Z0",  "f8",  "nT",  decimal_digits = 2, plot = True )        
    baseline_components += ["Z0"]
    
    add_result_comp(True, "M-offset",     "f8", "nT", plot = (calc_params is not None) and  calc_params.optimize_M_offset, decimal_digits = 3 )
    add_result_comp(True, "M-scale",      "f8", "",   plot = (calc_params is not None) and calc_params.optimize_M_scale,  decimal_digits = 5 )
    add_result_comp(True, "epsilon",      "f8", "arcmin",  plot = (calc_params is not None) and calc_params.optimize_epsilon,  decimal_digits = 3 )
    add_result_comp(True, "delta",        "f8", "arcsec", plot = (calc_params is not None) and calc_params.optimize_delta,    decimal_digits = 1 )
    
    add_result_comp( True, "M-offset_D/I-diff", "f8", "nT",     decimal_digits = 2,  plot = True)    
    add_result_comp( True, "epsilon_D/I-diff",  "f8", "arcmin", decimal_digits = 2,  plot = True)    
    
    add_result_comp(True, "DIM-I",       "f8", '°',  decimal_digits = angle_utils.angle_expression_decimal_digits['deg'], optimized_param = (calc_params is not None) and (not calc_params.method in [ "LM-D", "LM-|"]), plot = True )
    add_result_comp(True, "DIM-D",       "f8", '°',  decimal_digits = angle_utils.angle_expression_decimal_digits['deg'], optimized_param = (calc_params is not None) and (not calc_params.method in [ "LM-I", "LM-|"]), plot = True )
    
    add_result_comp(True, "DIM-X",       "f8", "nT", decimal_digits = 2, plot = True)
    add_result_comp(True, "DIM-Y",       "f8", "nT", decimal_digits = 2, plot = True)
    add_result_comp(True, "DIM-Z",       "f8", "nT", decimal_digits = 2, plot = True)
    # add_result_comp(True, "DIM-H",       "f8", "nT", decimal_digits = 2, plot = True)
    add_result_comp(True, "F",           "f8", "nT", decimal_digits = 2, plot = True)
    
    
    
    empty_temp_ch = 1
    while len(temp_comps) < 2: 
        while f"T{empty_temp_ch}" in temp_comps:
            empty_temp_ch += 1
        temp_comps += [ f"T{empty_temp_ch}" ]
        
   
    add_result_comp(False, temp_comps[0],  "f4", "degC",  decimal_digits = 2, plot = True)
    add_result_comp(False, temp_comps[1],  "f4", "degC",  decimal_digits = 2, plot = True)
    
    add_result_comp(False, "azimark_WGS84",          "f8", "deg",    decimal_digits = angle_utils.angle_expression_decimal_digits["deg"],  plot = False)  
    # add_result_comp(False, "azimark_hc_theo",        "f8", "deg",    decimal_digits = angle_utils.angle_expression_decimal_digits["deg"],  plot = False)  
    add_result_comp(True,  "azimark_hc-FL-FR-diff",  "f8", "arcsec", decimal_digits = 1,  plot = True )
    add_result_comp(True,  "azimark_vc-FL-FR-diff",  "f8", "arcsec", decimal_digits = 1,  plot = True )
    add_result_comp(True,  "vc_offset_from_azimark", "f8", "deg",    decimal_digits = 5,  plot = False )
    
    add_result_comp(False,  "F-pillar-diff", "f8", "nT",    decimal_digits = 2,  plot = False )  # It is reported by MagPy
    
   
    #add_result_comp(True, True, "Horiz-diff",  "f8",  "nT", plot = True )
    #if gu.COMPNAME_INCL_ANGLE in variometer_mag_data.dtype.names:
    #    add_result_comp(True, True, "Incl-diff", "f8", '"', plot = True )  
    #if gu.COMPNAME_DECL_ANGLE in variometer_mag_data.dtype.names:
    #    add_result_comp(True, True, "Decl-diff", '"',  plot = True )  #

    #if 'T1' in variometer_mag_data.dtype.names:
    #    add_result_comp(True, True, "Variometer-T1", "f8", 'K',  plot = True )
    #    add_result_comp(True, True, "Variometer-T2", "f8", 'K',  plot = True )

    return result_array_dtypes, result_array_comps_params, concatenated_plot_y_components


# Calculates spot baseline values for a variometer recording using D/I results of an absolute measurement (processed DIM data) 

def calc_baseline_values(calc_approach, mag_obs_dataset, variometer_mag_data, F_abs_pillar, DIM_magfield_Dpart, DIM_magfield_Ipart, common_calc_params):

    assert calc_approach in gu.BASELINE_CALCULATION_APPROACHES[1:]  # First possibility of "auto" is not accepted here
    assert isinstance(mag_obs_dataset, geomag_obs_dataset.mag_obs_dataset)

    # Transferring ADS's magfield for the variometer pillar

    DIM_magfield_at_variometer__D = copy.deepcopy(DIM_magfield_Dpart)
    DIM_magfield_at_variometer__I = copy.deepcopy(DIM_magfield_Ipart)
    if common_calc_params.dI_pillar_diff is not None:
        assert common_calc_params.dI_pillar_diff == 0, "Revise baseline calculation"
        DIM_magfield_at_variometer__D.I -= common_calc_params.dI_pillar_diff
        DIM_magfield_at_variometer__I.I -= common_calc_params.dI_pillar_diff
    if common_calc_params.dD_pillar_diff is not None:
        assert common_calc_params.dD_pillar_diff == 0, "Revise baseline calculation"
        DIM_magfield_at_variometer__D.D -= common_calc_params.dD_pillar_diff
        DIM_magfield_at_variometer__I.D -= common_calc_params.dD_pillar_diff
    if common_calc_params.F_pillar_diff_variometer is not None:
        assert common_calc_params.F_pillar_diff_variometer == 0, "Revise baseline calculation"
        DIM_magfield_at_variometer__D.F -= common_calc_params.F_pillar_diff_variometer
        DIM_magfield_at_variometer__I.F -= common_calc_params.F_pillar_diff_variometer

    mag_components_from_DIM = {}
    VAR_baseline_data = {}
    baseline_values   = {}

    HX_compname = mag_obs_dataset.usage_frame[0] 
    DY_compname = mag_obs_dataset.usage_frame[1] 


    if calc_approach == "Y0-nT":
        decbas = mag_obs_dataset.get_decbas()
        mag_components_from_DIM [HX_compname] =  DIM_magfield_at_variometer__I.calc_horiz_comp(decbas)
        # until 2025.12.15: mag_components_from_DIM [DY_compname] =  DIM_magfield_at_variometer__I.calc_horiz_comp(decbas + np.pi/2.0)
        mag_components_from_DIM [DY_compname] =  DIM_magfield_at_variometer__D.calc_horiz_comp(decbas + np.pi/2.0)
        VAR_baseline_data[HX_compname] = np.mean( variometer_mag_data[HX_compname][DIM_magfield_at_variometer__I.indeces] )
        # until 2025.12.15: VAR_baseline_data[DY_compname] = np.mean( variometer_mag_data[DY_compname][DIM_magfield_at_variometer__I.indeces] )
        VAR_baseline_data[DY_compname] = np.mean( variometer_mag_data[DY_compname][DIM_magfield_at_variometer__D.indeces] )

        baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR] = mag_components_from_DIM [HX_compname] - VAR_baseline_data[HX_compname]
        baseline_values[gu.COMPNAME__BASELINE__DY_SENSOR] = mag_components_from_DIM [DY_compname] - VAR_baseline_data[DY_compname]

    else:

        variometer_E__Incl_series  = np.mean( variometer_mag_data[DY_compname][DIM_magfield_at_variometer__I.indeces] )- 0*mag_obs_dataset.get_D0_nT_applied_at_file_reading()
        variometer_E__Decl_series  = np.mean( variometer_mag_data[DY_compname][DIM_magfield_at_variometer__D.indeces] )- 0*mag_obs_dataset.get_D0_nT_applied_at_file_reading()

        N_var = np.mean( variometer_mag_data[HX_compname][DIM_magfield_at_variometer__I.indeces] ) # - *mag_obs_dataset.get_H0_nT_applied_at_file_reading()

        F_Inkl_Avg_at_variometer = np.mean(F_abs_pillar[DIM_magfield_at_variometer__I.indeces]) - common_calc_params.F_pillar_diff_variometer

        N_b = baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR] =   np.sqrt( np.square(F_Inkl_Avg_at_variometer * np.cos(DIM_magfield_at_variometer__I.I) ) -  np.square( variometer_E__Incl_series )) - N_var

        mag_components_from_DIM['D'] =  DIM_magfield_at_variometer__D.D
        # D_var = VAR_baseline_data['D']  = np.arcsin( VAR_baseline_data['E-Decl']/DIM_magfield_at_variometer__D.H()  )
        D_variometer = np.arctan2( variometer_E__Decl_series,  np.mean(variometer_mag_data[HX_compname][DIM_magfield_at_variometer__D.indeces])-mag_obs_dataset.get_D0_nT_applied_at_file_reading() +  N_b)

        baseline_values[gu.COMPNAME__BASELINE__D_ANGLE] = mag_components_from_DIM['D'] -  D_variometer

        # Elvileg ez is jo, de rossz eredmény jön ki ezzel: mag_components_from_DIM ['Z']   = np.mean(F_abs_pillar * np.sin(DIM_magfield_at_variometer__I.I)) # QQQ DIM_magfield_at_variometer__I.Z()
        # mag_components_from_DIM ['Z'] = DIM_magfield_at_variometer__I.Z()
        # assert abs(old_Z-mag_components_from_DIM ['Z']) < 0.1

    VAR_baseline_data['Z']    = np.mean(variometer_mag_data['Z'][ DIM_magfield_at_variometer__I.indeces ])
    baseline_values[gu.COMPNAME__BASELINE__Z_SENSOR]      = DIM_magfield_at_variometer__I.Z() - VAR_baseline_data['Z']

    return baseline_values



# ---------------------------------------------------------------------------------------------------
# Conventional inversion of an absolute measurement program (processing window) using Levemberg-Marquardt algorithm.
# By optimization, it searches the D/I angles, and sensor errors (offset/scale/epsilon/delta).
# Absolute measurement data is in function parameters of 'i_mag_measured', 'i_hc_measured' / 'i_vc_measured'.
# On input, 'magfield' is the starting magnetic field vector (inclination/declination/F).
# On output, 'magfield'  contains the optimized results.
# The same stands for 'calc_var.instr': it is both the input and the output of the optimization.
# 'magfield' has to have a valid "I" value for the "known_inclination" sub-method (used to processdeclination-only measurement)
# ---------------------------------------------------------------------------------------------------
def mag_LM_inversions_for_DI_angles(i_mag_measured, mag_varcorrected, i_hc_measured_theo, i_vc_measured_theo, F,
                            magfield,
                            calc_var,
                            b_plot_best_result_residuals,
                            hc_tick_raw = None, b_print_new_random_step_results = True,
                            known_results = None,  # For simulation, it is a list of '[ magfield_class, instrument_param_class ]'
                            meas_method_str = "unknown",
                            variometer_data = None,
                            known_inclination = False, known_declination = False,
                            ):
    import lmfit

    assert isinstance(calc_var, spot_baselines_calculation_parameters_class)
    assert not np.any( np.isnan(i_mag_measured)  )

    dim_azimuth.set_default_azimuth_def(calc_var.true_angle_calculator)
    mag_sensor_index = calc_var.mag_sensor_index
    assert calc_var.mag_sensor_index>=0 and calc_var.mag_sensor_index<calc_var.instr.nof_mag_sensors()
    # b_optimize_scale = None # Global variable in this function
    # res_instr = copy.deepcopy(calc_var.instr) # Note: at the very end, new instrument params will be returned through 'calc_var.instr'
    # res_magfield = du.magfield_class()
    # calc_var_copy = copy.deepcopy(calc_var)  # Other name of this variable: "calculation_variation"

    model_params = None # We will use later assigned 'model_params' variable in 'process_results_and_format_solution_str()'


    n = len(i_mag_measured)

    assert  n == len(i_hc_measured_theo), "len(mag_measured) != len(hc_measured)"
    assert  n == len(i_vc_measured_theo), "len(mag_measured) != len(vc_measured)"

    mag_measured      = np.copy(i_mag_measured) + mag_varcorrected
    hc_measured_theo  = np.copy(i_hc_measured_theo)
    vc_measured_theo  = np.copy(i_vc_measured_theo)

    magfield.reference_coord_system = 'theo'
    # magfield.F = np.mean( F )

    # We can skip orto error calculation (for faster calculation) if fixed ortho error is 0.0
    # ToDo: check 'instr.ortho_error_arcsec == 0' expression agains numerical errors
    neglect_ortho_error = (not (calc_var.optimize_ortho_error))  and (calc_var.instr.ortho_error_arcsec == 0)
    assert neglect_ortho_error,  "DI modell neglects ortho error!"

    # hc_corr_functions_start = []
    # vc_corr_functions_start = []
    # for func in hc_corr_functions: hc_corr_functions_start = hc_corr_functions_start + copy.deepcopy(func)
    # for func in vc_corr_functions: vc_corr_functions_start = vc_corr_functions_start + copy.deepcopy(func)

    # if not AE1_hc_cos_function.ampl_optimize:
    #     instr.ae1_hc_ampl_arcmin = rad_to_arcmin( AE1_hc_cos_function.ampl )
    # if not AE1_vc_cos_function.ampl_optimize:
    #     instr.ae1_vc_ampl_arcmin = rad_to_arcmin( AE1_vc_cos_function.ampl )
    # if not AE1_hc_cos_function.phase_optimize:
    #     instr.ae1_hc_phase_deg = AE1_hc_cos_function.phase__as_deg()
    # if not AE1_vc_cos_function.phase_optimize:
    #     instr.ae1_vc_phase_deg = AE1_vc_cos_function.phase__as_deg()

    # ae1_hc_start.ampl = arcmin_to_rad ( instr.ae1_hc_ampl_arcmin )
    # assert ae1_hc_start.ampl >= arcsec_to_rad ( ae1_hc_start.ampl_arcsec_min)
    # assert ae1_hc_start.ampl <= arcsec_to_rad ( ae1_hc_start.ampl_arcsec_max)

    # ae1_vc_start.ampl = arcmin_to_rad ( instr.ae1_vc_ampl_arcmin )
    # assert ae1_vc_start.ampl >= arcsec_to_rad ( ae1_vc_start.ampl_arcsec_min)
    # assert ae1_vc_start.ampl <= arcsec_to_rad ( ae1_vc_start.ampl_arcsec_max)

    # ae1_hc_start.phase_ = instr.ae1_hc_phase_
    # ae1_vc_start.phase_ = instr.ae1_vc_phase_

    # ae2_start = copy.deepcopy(plate_common_error_cos_function)
    # if ae2_start.ampl_optimize:
    #     # Using instrument's plate error
    #     ae2_start.ampl = arcmin_to_rad(instr.ae2_ampl_arcmin)
    # if ae2_start.omega_optimize:
    #     ae2_start.omega = instr.ae2_omega
    # assert ae2_start.ampl >= arcsec_to_rad ( ae2_start.ampl_arcsec_min )
    # assert ae2_start.ampl <= arcsec_to_rad ( ae2_start.ampl_arcsec_max )

    # if ae2_start.phase_optimize:
    #     ae2_start.phase_ =  instr.ae2_phase_

    # if ae3_cos_function == None:
    #     ae3_cos_function_start = copy.deepcopy(angle_corr.no_120anglecorr)
    # else:
    #     ae3_cos_function_start = copy.deepcopy(ae3_cos_function)

    model = lmfit.Model(du.fullmagfit_diminstrument_based_DI_model, independent_vars = [ 'X', 'model_extra_params' ]  ) # varcorr_method=calc_var.varcorr_method, varcorr_data=varcorr_data

    # Initial values: using received  'magfield' and 'instr' values

    mag_sensor =  calc_var.instr.mag_sensors[mag_sensor_index]
    assert isinstance(mag_sensor, dim_instr.mag_sensor_model_class)

    model_params = lmfit.Parameters()
    calc_var.add_default_instr_model_params(model_params)

    model_params.add_many(
        lmfit.Parameter('I', value = magfield.I, min = magfield.I - math.radians(calc_var.max_incl_change_deg),    max = magfield.I + math.radians(calc_var.max_incl_change_deg), user_data = model_utils.conv_radian_to_deg, vary = not known_inclination ),
        lmfit.Parameter('D', value = magfield.D, min = magfield.D - math.radians(calc_var.max_decl_change_deg),    max = magfield.D + math.radians(calc_var.max_decl_change_deg), user_data = model_utils.conv_radian_to_deg, vary = not known_declination ),
    )

    # Moved into 'model_extra_params': model.set_param_hint('neglect_ortho_error',      value = neglect_ortho_error, vary = False)

    # model.set_param_hint('ae1_hc_ampl',    value = ae1_hc_start.ampl,  min = arcsec_to_rad(ae1_hc_start.ampl_arcsec_min) ,   max = arcsec_to_rad(ae1_hc_start.ampl_arcsec_max), vary = ae1_hc_start.ampl_optimize)
    # model.set_param_hint('ae1_vc_ampl',    value = ae1_vc_start.ampl,  min = arcsec_to_rad(ae1_vc_start.ampl_arcsec_min) ,   max = arcsec_to_rad(ae1_vc_start.ampl_arcsec_max), vary = ae1_vc_start.ampl_optimize)

    # # ToDo: use 'ae1_hc_start.phase_optimize' instead of 'calc_var.optimize_excentr_error'
    # model.set_param_hint('ae1_hc_phase_',   value = ae1_hc_start.phase_,   min = instr.ae1_hc_phase_min,   max = instr.ae1_hc_phase_max, vary = ae1_hc_start.phase_optimize)
    # model.set_param_hint('ae1_vc_phase_',   value = ae1_vc_start.phase_,   min = instr.ae1_vc_phase_min,   max = instr.ae1_hc_phase_max, vary = ae1_vc_start.phase_optimize)

    # model.set_param_hint('ae2_ampl',   value = ae2_start.ampl,   min = arcsec_to_rad(ae2_start.ampl_arcsec_min) ,  max = arcsec_to_rad(ae2_start.ampl_arcsec_max),  vary = ae2_start.ampl_optimize)
    # model.set_param_hint('ae2_phase_', value = ae2_start.phase_, min = MIN_PHASE_,                                                max = MAX_PHASE_,                                               vary = ae2_start.phase_optimize)
    # model.set_param_hint('ae2_omega',  value = ae2_start.omega,  min = instr.ae2_omega_min,            max = instr.ae2_omega_max,           vary = ae2_start.omega_optimize)

    # model.set_param_hint('ae3_ampl',   value = ae3_cos_function_start.ampl,   min = arcsec_to_rad(ae3_cos_function_start.ampl_arcsec_min) , max = arcsec_to_rad(ae3_cos_function_start.ampl_arcsec_max),  vary = ae3_cos_function_start.ampl_optimize)
    # model.set_param_hint('ae3_phase_', value = ae3_cos_function_start.phase_, min = MIN_PHASE_,                                             max = MAX_PHASE_,                                             vary = ae3_cos_function_start.phase_optimize)
    # model.set_param_hint('ae3_omega',  value = ae3_cos_function_start.omega,  min = ae3_cos_function_start.omega_min,                       max = ae3_cos_function_start.omega_max,                       vary = ae3_cos_function_start.omega_optimize)
    # model.set_param_hint('consider_ae3' , value = ae3_cos_function_start.is_active(), vary = False )



    if (meas_method_str=='NM'): #  (meas_method_str[0:8] == "nullmeth"):
        # In a "standard" null method measurement, we have eight measurement parts within a null method measurement, plus the scale check
        # We use this info in the plot:
        meas_point_groups = du.NM_separate_points_groups(i_hc_measured_theo, i_vc_measured_theo, calc_var.instr)
    else:
        meas_point_groups = None

    # Running optimization:
    modelfunc_X =   [ hc_measured_theo, vc_measured_theo, F ]
    safety_check__angles_does_not_change = hc_measured_theo[0]

    if calc_var.varcorr_method=="LM/ID-var":
        # varcorr_method_index = 1
        assert variometer_data is not None
        if (gu.COMPNAME_INCL_ANGLE in variometer_data.dtype.names) and (gu.COMPNAME_DECL_ANGLE in variometer_data.dtype.names):
            I_varcorr = variometer_data[gu.COMPNAME_INCL_ANGLE]
            D_varcorr = variometer_data[gu.COMPNAME_DECL_ANGLE]
        else:
            D_varcorr, I_varcorr = gu.calc_D_I_angles(variometer_data['H'],variometer_data['D'],variometer_data['Z'] )

        modelfunc_X += [ np.mean(I_varcorr) - I_varcorr , np.mean(D_varcorr)  - D_varcorr    ]

    model_extra_params = {'varcorr_method':calc_var.varcorr_method, 'variometer_data':variometer_data, 'F_data': F,  "neglect_ortho_error": neglect_ortho_error, "M_DIM":None }
    
    if (calc_var.test_mode is not None) and calc_var.test_mode:
        model_extra_params["M_DIM"] = mag_measured
        target_values = np.zeros(len(mag_measured))
    else:
        target_values =  mag_measured

    best_resid_std, best_results, params_at_the_bound_names, bestres_random_try_no, bestres_nof_removed_outliers, best_res_residuals = model_utils.lmfit_with_randominitialvalues_execute(calc_var.random_steps, model, model_params, target_values, modelfunc_X, reference_params = None, reference_params_descr = "", reference_params_notused_indeces = [],  b_plot_best_result_residuals = b_plot_best_result_residuals, b_print_new_random_step_results = b_print_new_random_step_results, model_extra_params =  model_extra_params, meas_point_groups=meas_point_groups)

    assert safety_check__angles_does_not_change == hc_measured_theo[0]

    # Best results' are now stored in 'res_instr_params' and 'res_magfield'



    #   best_result_str =  process_results_and_format_solution_str( model_params,   best_resid_std, best_results, "oneline" )
    best_result_str =  model_utils.gen_get_solution_str_callback(model_params, best_resid_std, best_results)
    best_result_str =  " {:02d} {:4d} {:d} |  ".format(n-bestres_nof_removed_outliers, bestres_random_try_no, bestres_nof_removed_outliers) # + best_result_str

    assert model_params['I'].value == best_results[9]
    assert model_params['D'].value == best_results[10]

    # This is a output object:
    magfield.F = np.mean( F )
    magfield.D = model_params['D'].value
    magfield.I = model_params['I'].value

    # Copying results into 'instr' that is sent back to the caller
    # 2023.04.29: we copy the offset/scale/... parameters even if they are not optimized, as the fixed starting offset/scale/.... parameters can be other than the values stored in instr

    # assert abs( instr.mag_sensors[mag_sensor_index].offset - res_instr.mag_sensors[mag_sensor_index].offset )  < 1e-5
    # assert abs( instr.mag_sensors[mag_sensor_index].scale  - res_instr.mag_sensors[mag_sensor_index].scale) < 1e-5

    calc_var.apply_instr_model_params(model_params)

    # if True or calc_var.optimize_M_offset:    instr.mag_sensors[mag_sensor_index].offset           = res_instr.mag_sensors[mag_sensor_index].offset
    # if True or calc_var.optimize_M_scale:     instr.mag_sensors[mag_sensor_index].scale            = res_instr.mag_sensors[mag_sensor_index].scale
    # if True or calc_var.optimize_epsilon:     instr.mag_sensors[mag_sensor_index].epsilon_arcmin   = res_instr.mag_sensors[mag_sensor_index].epsilon_arcmin
    # if True or calc_var.optimize_delta:       instr.mag_sensors[mag_sensor_index].delta_arcsec     = res_instr.mag_sensors[mag_sensor_index].delta_arcsec
    # if True or calc_var.optimize_ortho_error: instr.ortho_error_arcsec = res_instr.ortho_error_arcsec

    if calc_var.optimize_hc_eccentr_err:
        assert False,"Implement these in 'apply_instr_model_params:"
        # instr.hc.eccentr_error_phase = res_instr.hc.eccentr_error_phase
        # instr.hc.eccentr_error_ampl = res_instr.hc.eccentr_error_ampl
    if calc_var.optimize_vc_eccentr_err:
        assert False,"Implement these in 'apply_instr_model_params:"
        # instr.vc.eccentr_error_phase = res_instr.vc.eccentr_error_phase
        # instr.vc.eccentr_error_ampl = res_instr.vc.eccentr_error_ampl

    # if True or calc_var.optimize_ortho_error: instr.ortho_error_arcsec = res_instr.ortho_error_arcsec

    return best_result_str, params_at_the_bound_names, bestres_nof_removed_outliers, best_resid_std, best_res_residuals


# ---------------------------------------------------------------------------------------------------
# This type of inversion estimates average baseline values + sensor-parameter values + optionally decbas .
# Suitable for per-processing-window and for and "All-cycles-in-one" optimiziations.
# Baselines and decbas are calculated for 'variometer_data_def'
# Three submethods:
#  - "A0B0":  only two baselines are optimized, the third is determined from X_abs^2 + Y_abs^2 + Z_abs^2 = F^2 relation
#  - "A0B0Z0": all three baselines are optimized
#  - "H0D0": two baselines are optimized. Similar approach than null-method solver uses for HDZ variometer
# For input data, there are two main options:
#  - input data has reference position data, so eccentricity errors can be estiamted 
#  - input data has no reference positions data, so eccentricity can not be estimated
#
# Megj.: előfordulhatnak  a valóságossal ekvivalens nem reális, de jó megoldások.
#   Az eredő F/I/D magfield ilyenkor is jól becsült lesz, és a később, a hívó rutinban, az F/I/D-ből kiszámolt A0/B0/Z0 már jó értékű lesz.
# ---------------------------------------------------------------------------------------------------
def mag_LM_inversions_for_baseline(baseline_calc_approach, has_refpos_data, sub_method, input_data, calc_var_index, calc_var, mag_obs_dataset, DIM_instr, b_plot_best_result_residuals, variometer_mag_data, b_print_new_random_step_results = True, b_plot__F_diff_vs_posID = True, show_results_table = True):
    import lmfit
    assert sub_method in ["A0B0", "A0B0Z0", "H0D0" ]

    assert isinstance(calc_var, spot_baselines_calculation_parameters_class)
    assert isinstance(mag_obs_dataset, geomag_obs_dataset.mag_obs_dataset)
    calcvar_index_str = calc_var.id
    if has_refpos_data:
        # Originally this is the case when ALL-IN-ONE optimization is run
        # 'ads_mag_data_with_refpos' is a structured array, returned by 'deprecated_get_all_records_with_refpos_data()'
        ads_mag_data_with_refpos = input_data # It is a structured array
        mag_measured     = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__MAG]
        assert not np.any( np.isnan(mag_measured)  )

        # hc_measured_theo = calc_var.correct_raw_angles( "hc",  dim_instr.tick_to_rad(ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__HC_TICK], False) )
        # vc_measured_theo = calc_var.correct_raw_angles( "vc",  dim_instr.tick_to_rad(ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__VC_TICK], True, calc_var.instr.vc0) )

        # hc_MIRA_SU_theo_measured = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__REFPOS_HC_THEO_SU]
        # hc_MIRA_SD_theo_measured = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__REFPOS_HC_THEO_SD]
        # hc_MIRA_SU_WGS84 = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__REFPOS_HC_WGS84_SU]
        # hc_MIRA_SD_WGS84 = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__REFPOS_HC_WGS84_SD]

        timelabels = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__TIMESTAMP]

        recalc_NED_coordinates = False #
        if not recalc_NED_coordinates:
            precalc_hc_NED = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__PRECALC_HC_NED_+calcvar_index_str] # dim_azimuth.h_theo_to_ned(hc_measured_theo, hc_MIRA_SU_theo_measured, hc_MIRA_SD_theo_measured, hc_MIRA_SU_WGS84, hc_MIRA_SD_WGS84 )
            precalc_vc_NED = ads_mag_data_with_refpos[dim_measurement_loader.COMPNAME__PRECALC_VC_NED_+calcvar_index_str] #  dim_azimuth.v_theo_to_ned(vc_measured_theo)
    else:
        timelabels       = input_data[0]
        mag_measured     = input_data[1]
        # precalc_hc_theo  = input_data[2]  # Or: dim_azimuth.default_azimuth_def.h_theo_to_ned(hc_measured_theo)
        # precalc_vc_theo  = input_data[3]
        recalc_NED_coordinates = False
        precalc_hc_NED =  input_data[4]
        precalc_vc_NED =  input_data[5]
        F_data         =  input_data[6]

    # assert isinstance(variometer_data_def, du.mag_data_set_def)
    
    assert baseline_calc_approach == "Y0-nT", "Not supported baseline-type for LM calculation: {}".format(baseline_calc_approach)
    
    calc_frame = mag_obs_dataset.usage_frame

    assert calc_frame in [ "HDZ", "XYZ" ]

    # Not required, as 'ads_mag_data_with_refpos' has the preprocessed MIRA data:  dim_azimuth.set_default_azimuth_def(calc_var.azimuth_def)
    mag_sensor_index = calc_var.mag_sensor_index

    # res_instr = copy.deepcopy(instr) # Note: at the very end, new instrument params will be returned through 'instr'
    # res_magfield = du.magfield_class()
    # calc_var = copy.deepcopy(calc_var)  # Other name of this variable: "calculation_variation"

    # calc_F = mag_obs_dataset.variometer_data.calc_F(variometer_mag_data)
    # plot_utils.plot_datetime_y( timelabels,  [calc_F - F_data ], "delta-F before", fmt = F_DIFF_FMT ,color = F_DIFF_COLOR )

    model_params = None # We will use later assigned 'model_params' variable in 'process_results_and_format_solution_str()'

    def ads_variometer_baseline_diff_model(X, baseline_A0, baseline_B0, baseline_Z0,
                                              decbas, delta_F,
                                              M_scale, M_offset, epsilon, delta,
                                              hc_ecc_phase, vc_ecc_phase,  hc_ecc_ampl,  vc_ecc_ampl,
                                              ortho_error, model_extra_params = []):
        var_A = X[0]  # H or X
        var_B = X[1]  # D or Y
        var_Z = X[2]  # Always Z
        var_F = X[3]

        # hc_measured_theo = X[4]
        # vc_measured_theo = X[5]

        # F should be a time series:
        # assert len(var_F) == len(hc_measured_theo)

        # nonlocal recalc_NED_coordinates
        # if recalc_NED_coordinates:
        #     hc_MIRA_SU_theo_measured = X[6]
        #     hc_MIRA_SD_theo_measured = X[7]
        #     hc_MIRA_SU_WGS84 = X[8]
        #     hc_MIRA_SD_WGS84 = X[9]
        #     if np.isnan(hc_ecc_ampl):
        #         hc_theo = hc_measured_theo
        #     else:
        #         hc_theo = hc_measured_theo + hc_ecc_ampl * np.cos( hc_measured_theo + hc_ecc_phase )
        #     hc_MIRA_SU = hc_MIRA_SU_theo_measured + hc_ecc_ampl * np.cos( hc_MIRA_SU_theo_measured + hc_ecc_phase )
        #     hc_MIRA_SD = hc_MIRA_SD_theo_measured + hc_ecc_ampl * np.cos( hc_MIRA_SD_theo_measured + hc_ecc_phase )
        #     hpa_n = dim_azimuth.h_theo_to_ned(hc_theo, hc_MIRA_SU,hc_MIRA_SD, hc_MIRA_SU_WGS84, hc_MIRA_SD_WGS84 )
        #     if not np.isnan(vc_ecc_ampl):
        #         vc_theo = vc_measured_theo + vc_ecc_ampl * np.cos( vc_measured_theo + vc_ecc_phase )
        #     else:
        #         vc_theo = precalc_vc_theo
        # else:
        #     nonlocal calc_var
        #     assert (not calc_var.optimize_hc_eccentr_err) and  (not calc_var.optimize_vc_eccentr_err)
        #     hpa_n = precalc_hc_theo
        #     vpa_n = precalc_vc_theo

        # if not np.isnan(vc_ecc_ampl):
        #     vc_theo = vc_measured_theo + vc_ecc_ampl * np.cos( vc_measured_theo + vc_ecc_phase )
        # else:
        #     vc_theo = vpa_n = precalc_vc_theo
        # vpa_n = dim_azimuth.v_theo_to_ned(vc_theo)

        # print("F:{:.1f} hc={:.1f}° vc={:.1f}°  I={:.2f}° D={:.2f}°  epsilon:{:.1f}' delta:{:.1f}'  Moffs={:.1f} nT scale={:.1f} nT   ".format(F, math.degrees(hc_measured[0]), math.degrees(vc_measured[0]), math.degrees(I), math.degrees(D), rad_to_arcsec(epsilon), rad_to_arcsec(delta), Moffset, scale ))
        baselinecorrected_fluxgate_A = var_A + baseline_A0
        baselinecorrected_fluxgate_B = var_B + baseline_B0
        nonlocal sub_method
        if sub_method=="A0B0Z0":
            baselinecorrected_fluxgate_Z = var_Z + baseline_Z0  # Z+Z0
        else:
            baselinecorrected_fluxgate_Z = np.mean( np.sqrt( np.square(var_F + delta_F) - np.square(baselinecorrected_fluxgate_A) - np.square(baselinecorrected_fluxgate_B)   ) )
            nonlocal model_params
            assert not model_params['baseline_Z0'].vary
            model_params['baseline_Z0'].value = baselinecorrected_fluxgate_Z - np.mean(var_Z)


        D,I = gu.calc_D_I_angles(decbas, baselinecorrected_fluxgate_A,baselinecorrected_fluxgate_B,baselinecorrected_fluxgate_Z)

        assert not np.isnan(I).any()
        nonlocal precalc_hc_NED, precalc_vc_NED
        hpa_n = precalc_hc_NED
        vpa_n = precalc_vc_NED
        return amdp.dim_instrument_model(var_F + delta_F, vpa_n, hpa_n, I, D, epsilon, delta, M_offset, M_scale)

    n = len(mag_measured)
    assert  n == len(precalc_hc_NED), "len(mag_measured) != len(hc_measured)"
    assert  n == len(precalc_vc_NED), "len(mag_measured) != len(vc_measured)"

     # We can skip orto error calculation (for faster calculation) if fixed ortho error is 0.0
    neglect_ortho_error = True # (not (calc_var.optimize_ortho_error))  and (calc_var.instr.ortho_error_arcsec == 0)
    assert neglect_ortho_error,  "DI modell neglects ortho error!"

    model = lmfit.Model(ads_variometer_baseline_diff_model, independent_vars = [ 'X', 'model_extra_params' ] )

    mag_sensor =  DIM_instr.mag_sensors[mag_sensor_index]
    assert isinstance(mag_sensor, dim_instr.mag_sensor_model_class)

    model_params = lmfit.Parameters()

    # if variometer_data_def.frame == "XYZ":
    #     starting_decbas = 0.0
    # elif variometer_data_def.frame == "HDZ":
    #     starting_decbas = math.radians(measurement_sites.DECL_START_DEG) # average declination
    # else:
    #     assert False

    starting_decbas = mag_obs_dataset.get_decbas()
    max_decbas_change = np.radians(2.5)

    model_params.add_many(
        lmfit.Parameter('baseline_A0',   value = calc_var.opt_start["BL-A0"],             min =  calc_var.opt_range["BL-A0"][0],        max = calc_var.opt_range["BL-A0"][1],          vary = True, user_data = ("nT",)),
        lmfit.Parameter('baseline_B0',   value = calc_var.opt_start["BL-B0"],              min =  calc_var.opt_range["BL-B0"][0],        max = calc_var.opt_range["BL-B0"][1],          vary = True, user_data = ("nT",)),
        lmfit.Parameter('baseline_Z0',   value = calc_var.opt_start["BL-Z0"],              min =  calc_var.opt_range["BL-Z0"][0],        max = calc_var.opt_range["BL-Z0"][1],          vary = sub_method=="A0B0Z0",  user_data = ("nT",)),
        lmfit.Parameter('decbas',        value = starting_decbas,   min = starting_decbas-max_decbas_change,   max = starting_decbas+max_decbas_change,    vary = calc_var.optimize_decbas,  user_data = model_utils.conv_radian_to_deg ),
        lmfit.Parameter('delta_F',       value = 0.0,               min = -20,                                 max = +20,        vary = False,  user_data = ("nT",))
     )

    calc_var.add_default_instr_model_params(model_params, DIM_instr)

    if  (calc_frame == "XYZ"):
        modelfunc_X =  [ variometer_mag_data['X'], variometer_mag_data['Y'], variometer_mag_data['Z'], F_data ]
    elif (calc_frame == "HDZ"):
        modelfunc_X = [ variometer_mag_data['H'], variometer_mag_data['D'], variometer_mag_data['Z'], F_data ]
    else:
        assert False

    # modelfunc_X +=   [ hc_measured_theo, vc_measured_theo ]
    # if has_refpos_data:  modelfunc_X +=  [ hc_MIRA_SU_theo_measured, hc_MIRA_SD_theo_measured, hc_MIRA_SU_WGS84, hc_MIRA_SD_WGS84 ]

    # Running optimization:
    # safety_check__angles_does_not_change = hc_measured_theo[0]
    model_extra_params = []

    best_resid_std, best_results, params_at_the_bound_names, bestres_random_try_no, bestres_nof_removed_outliers, best_res_residuals = \
        model_utils.lmfit_with_randominitialvalues_execute(calc_var.random_steps, model, model_params, mag_measured, modelfunc_X, reference_params = None, reference_params_descr = "", reference_params_notused_indeces = [], get_solution_str_callback = None, b_plot_best_result_residuals = b_plot_best_result_residuals, b_print_new_random_step_results = b_print_new_random_step_results, model_extra_params =  model_extra_params, timelabels = timelabels, show_results_table =  show_results_table)

    # assert safety_check__angles_does_not_change == hc_measured_theo[0]

    best_result_str =  model_utils.gen_get_solution_str_callback(model_params, best_resid_std, best_results)
    best_result_str =  " {:02d} {:4d} {:d} |  ".format(n-bestres_nof_removed_outliers, bestres_random_try_no, bestres_nof_removed_outliers) + best_result_str
    # print("Best result: {}".format(best_result_str) )

    calc_var.apply_instr_model_params(model_params, DIM_instr) # Updating 'calc_var.instr' to the optimized parameters

    model_utils.model_params_set_from_results(model_params, best_results)

    decbas                       = model_params['decbas'].value
    print("decbas= {:.4f} ° ".format(math.degrees(decbas)) )
    F_corr                       = F_data + model_params['delta_F'].value
    baseline_A0                  = model_params['baseline_A0'].value
    baseline_B0                  = model_params['baseline_B0'].value
    baselinecorrected_fluxgate_A = variometer_mag_data[ calc_frame[0] ] + baseline_A0  # H or X
    baselinecorrected_fluxgate_B = variometer_mag_data[ calc_frame[1] ] + baseline_B0  # D or Y
    if sub_method == "A0B0Z0":
        baseline_Z0                  = model_params['baseline_Z0'].value
    else:
        # 'model_params['baseline_Z0']' should contain the calculated Z0, however, for self-check, we recalcualte it:
        baseline_Z0   =  np.mean( np.sqrt( np.square(F_corr) - np.square(baselinecorrected_fluxgate_A)  - np.square(baselinecorrected_fluxgate_B) ) - variometer_mag_data['Z'] )
        if abs(model_params['baseline_Z0'].value - baseline_Z0 ) > 0.2: cu.warning( "Z0 self-check error: {} - {}".format(model_params['baseline_Z0'].value, baseline_Z0))
    baselinecorrected_fluxgate_Z = variometer_mag_data['Z'] + baseline_Z0

    D, I = gu.calc_D_I_angles(decbas, baselinecorrected_fluxgate_A, baselinecorrected_fluxgate_B, baselinecorrected_fluxgate_Z )

    calc_F = np.sqrt( np.square(baselinecorrected_fluxgate_A)  + np.square(baselinecorrected_fluxgate_B) + np.square(baselinecorrected_fluxgate_Z)  )

    if (b_plot__F_diff_vs_posID):
        F_diff = calc_F - F_corr
        plot_utils.plot_datetime_y( timelabels ,  [ F_diff  ], "delta-F after",fmt = F_DIFF_FMT , color = F_DIFF_COLOR )
        plot_utils.plot_xy( input_data[dim_measurement_loader.COMPNAME__POS_ID] , F_diff, "pos-ID", "delta-F after", symbol = F_DIFF_FMT , color = F_DIFF_COLOR )


    # Copying results into 'instr' that is sent back to the caller
    # 2023.04.29: we copy the offset/scale/... parameters even if they are not optimized, as the fixed starting offset/scale/.... parameters can be other than the values stored in instr

    return best_result_str, params_at_the_bound_names,  bestres_nof_removed_outliers, best_resid_std, best_res_residuals, [I,D,F_corr,decbas,baseline_A0,baseline_B0,baseline_Z0]


# ---------------------------------------------------------------------------------------------------
# Inversion of an absolute measurement (processing-window) using Gauss-Newton algorithm.
# ---------------------------------------------------------------------------------------------------
# 2025.03.22 előtt ezt kapta paraméterul: 'instr, mag_sensor_index, ' de elég neki a 'mag_sensor' mint paraméter!
def mag_inverse_method_with_GN(dim_mag_measured, dim_mag_varcorr, hc_NED, vc_NED, F, DIM_instr, calc_var, magfield, b_plot_best_result_residuals, variometer_mag_data, mag_obs_dataset):

    # assert DIM_instr.ortho_error_arcsec == 0.0
    assert (not calc_var.optimize_ortho_error) and (not calc_var.optimize_hc_eccentr_err) and (not calc_var.optimize_vc_eccentr_err)
    assert calc_var.optimize_M_scale and calc_var.optimize_M_offset and calc_var.optimize_epsilon and calc_var.optimize_delta

    mag_sensor = DIM_instr.mag_sensors[calc_var.mag_sensor_index]

    assert len(F) == len(hc_NED)

    MAX_ITER_COUNTER = 20
    iter_count = 0
    while True:

        regressors = np.vstack( (F, vc_NED, hc_NED) ).T
        parameters = (magfield.I, magfield.D,  mag_sensor.epsilon(), mag_sensor.delta(), mag_sensor.offset, mag_sensor.scale )

        parameters, stat, message, residuals = amdp.gauss_newton(dim_mag_measured + dim_mag_varcorr, regressors, parameters)
        magfield.I = parameters[0]
        magfield.D = parameters[1]
        magfield.F = np.mean(F)

        mag_sensor.epsilon_arcmin = angle_utils.rad_to_arcmin( parameters[2] )
        mag_sensor.delta_arcsec   = angle_utils.rad_to_arcsec( parameters[3] )

        mag_sensor.offset =  parameters[4]
        mag_sensor.scale =   parameters[5]

        if magfield.I > np.pi: magfield.I = magfield.I - mu.pi_x_2  # I ranges from -180.0 deg to +180.0 deg

        if iter_count == MAX_ITER_COUNTER: break
        # Recalculating correction values with the new results:
        baseline_values = calc_baseline_values(calc_var.baseline_calc_approach,  mag_obs_dataset, variometer_mag_data, F, magfield, magfield)
        prelim_F,prelim_I,prelim_D = gu.variometer_to_abs_data(calc_var.baseline_calc_approach, mag_obs_dataset, variometer_mag_data, baseline_values)
        new_dim_mag_varcorr = du.calc_absmeas_variation_correction_values(calc_var.varcorr_method, hc_NED, vc_NED, mag_sensor, F, prelim_I, prelim_D )
        change_offset, change_stddev = np.mean(new_dim_mag_varcorr-dim_mag_varcorr) , np.std(new_dim_mag_varcorr-dim_mag_varcorr)
        if (iter_count>0) and  (abs(change_offset)<0.05) and (abs(change_stddev)<0.02): break
        varcorr_mag_sensor = copy.copy(mag_sensor)
        # print("============ {:.2f}  {:.2f} =============".format(change_offset,change_stddev))
        dim_mag_varcorr = new_dim_mag_varcorr
        iter_count += 1


    return residuals, varcorr_mag_sensor, new_dim_mag_varcorr


# Processes residual method processing window with the analytical method.
# Based on paper K. Lauridsen "Experiences with the DI-fluxgate Magnetometer Inclusive Theory of the 
# Instrument and Comparison with Other Methods" (1985)" with note of: "This version is including the 
# corrections made by Kring Lauridsen in the printed version - 2008, Jürgen Matzka, DMI"
# 
# Optimized instrument params are returned on 'instr', while D and I are returned directly.
# 'F' is the F variation at the absolute pillar
# It can process four possible measurement data
#  1) Conventional null method with 8+1(optional scale check) meausrement
#  2) 8+1(optional scale check)  residual measuremett
#  3) (8*x)+1 "scan" measurement, that can be treated as
#      3/A) residual measruement
#      3/B) pseudo null measurement (deprecated)
def residual_method_anal_solver(baseline_calc_approach, position_indeces, measured_mag, hc, vc, hc_NED, vc_NED,F,instr,res_mag_sensor,test_mode, variometer_data, mag_obs_dataset, common_calc_params, calc_params, meas_pos):

    assert mag_obs_dataset.usage_frame in ["HDZ", "XYZ"]
    
    mag = np.copy(measured_mag)
    
    H0 = np.nan
    D0 = np.nan
    Z0 = np.nan

    MAX_ITERATIONS = 10 # 2...3: seems to be generally ok, but in a few cases, too small! 
                        # Bigger value sometime help result quality, but maybe sometime introduce error?!

    debug_iteration_process = False
    if debug_iteration_process:
        iteration_progress_data = np.empty( (MAX_ITERATIONS), dtype=[ ('I_I',float), ('D_I',float), ('delta_arcsec',float), ('M_scale',float)]  )

    mode = ["0nT-scan","true-resid"][1]

    delta_prelim =  0.0
    if mode == "true-resid":
        # 'I__I_part_avg_prelim' / 'D__I_part_avg_prelim' must be known
        if (calc_params is None) :
            I__I_part_avg_prelim = meas_pos.I_typical
            D__I_part_avg_prelim = meas_pos.D_typical
        else:
            I__I_part_avg_prelim = calc_params.last_results['I']
            D__I_part_avg_prelim = calc_params.last_results['D']
            if 'mag_sensor' in calc_params.last_results.keys():
                delta_prelim  = calc_params.last_results['mag_sensor'].delta()
    elif mode == "0nT-scan":
        assert False, "Deprecated"
        NM_meas_part_results, nr_of_scan_points, has_scale_check = du.NM_calculate_zeropos_angles_and_positions(hc,vc,mag,instr,calc_params.true_angle_calculator)
        if NM_meas_part_results == None:
            # Points could not be separated into at least 8 parts:
            cu.warning("Could not separate measurement points into Null-method measurement sections", True)
            return None, None, None, None, None, None, None, None, None, None
    else:
        assert False, "Unsupported mode: \"{}\"".format(mode)

    ED_indeces = position_indeces['ER']  # In a standard DIM instrument, the "Down" sensor position corresponds to "Face-Right" direction
    WD_indeces = position_indeces['WR']
    EU_indeces = position_indeces['EL']
    WU_indeces = position_indeces['WL']
    NU_indeces = position_indeces['NL']
    SD_indeces = position_indeces['SR']
    ND_indeces = position_indeces['NR']
    SU_indeces = position_indeces['SL']

    Decl_series_indeces = ED_indeces + WD_indeces + EU_indeces + WU_indeces
    Incl_series_indeces = NU_indeces + SD_indeces + ND_indeces + SU_indeces
    SC_indeces = position_indeces[du.SCALECHECK_POS_ID]
    
    if not calc_params.allow_unequal_measurement: 
        if len(Decl_series_indeces) != len(Incl_series_indeces): return "Unbalanced measurement", None
        
    assert len(Decl_series_indeces) % 4 == 0
    assert len(Decl_series_indeces) >= 4

    # theta<n> = vc angle, where 0 deg is horizontal SD
    # alpha<n> = hc_NED angle

    M_new_scale = 1.0
    b_estimate_M_scale = calc_params.optimize_M_scale and np.any( abs(mag) > 0.5 )    # Some manual measurements do not have scale value 
    
    if mag_obs_dataset is not None:
        H_comp_name  = mag_obs_dataset.usage_frame[0] # variometer_data.dtype.names[0] (not working for DIDD with comp order of F,D+,...)
        D_comp_name  = mag_obs_dataset.usage_frame[1] # variometer_data_def.instr_frame[1]


    for iteration_count in range(MAX_ITERATIONS):

        # --------------------------------------
        # --- Processing of inclination part ---
        # --------------------------------------

        F__I_part_avg = np.mean( F[ Incl_series_indeces ] )

        S1 = np.mean( mag[NU_indeces] )
        S2 = np.mean( mag[SD_indeces] )
        S3 = np.mean( mag[ND_indeces] )
        S4 = np.mean( mag[SU_indeces] )

                    
        MAX_ANGLE_ERROR = np.radians(35) if iteration_count<MAX_ITERATIONS/2 else np.radians(15)  

        
        try:  # 'select_equivalent_angle' drops an Exception in case of angle problem
            # alpha_xx: horizontal angles relative to the magnetic meridian:
            # 'np.mean()' should be outside of 'select_equivalent_angle() function!
            alpha1 = np.mean( angle_utils.select_equivalent_angle( hc_NED[NU_indeces]  - D__I_part_avg_prelim, 0.0,   MAX_ANGLE_ERROR, False ) )
            alpha2 = np.mean( angle_utils.select_equivalent_angle( hc_NED[SD_indeces]  - D__I_part_avg_prelim, 0.0,   MAX_ANGLE_ERROR, False ) )
            alpha3 = np.mean( angle_utils.select_equivalent_angle( hc_NED[ND_indeces]  - D__I_part_avg_prelim, np.pi, MAX_ANGLE_ERROR, False ) )
            alpha4 = np.mean( angle_utils.select_equivalent_angle( hc_NED[SU_indeces]  - D__I_part_avg_prelim, np.pi, MAX_ANGLE_ERROR, False ) )
            
            if np.isnan(alpha1) or np.isnan(alpha2) or np.isnan(alpha3) or np.isnan(alpha4):
                # This typically happens if reference mark angle is wrong
                return "Magnetic meridian angle error", None
                
            theta1 = np.mean( angle_utils.select_equivalent_angle( vc_NED[NU_indeces], + I__I_part_avg_prelim - np.pi/2,    MAX_ANGLE_ERROR) )
            theta2 = np.mean( angle_utils.select_equivalent_angle( vc_NED[SD_indeces], + I__I_part_avg_prelim + np.pi/2,    MAX_ANGLE_ERROR) )
            theta3 = np.mean( angle_utils.select_equivalent_angle( vc_NED[ND_indeces], 3*np.pi/2 - I__I_part_avg_prelim,   MAX_ANGLE_ERROR) )
            theta4 = np.mean( angle_utils.select_equivalent_angle( vc_NED[SU_indeces],   np.pi/2 - I__I_part_avg_prelim,   MAX_ANGLE_ERROR) )
        except ValueError:
            return "alpha or theta calculation error due to some angle problem", None    
            

        I__I_part_avg_ideal = (theta1 + theta2 - theta3 - theta4) / 4.0 + np.pi / 2.0 # amdp.Inclination(theta1,theta2,theta3,theta4)
        assert abs(angle_utils.angle_distance(I__I_part_avg_ideal, I__I_part_avg_prelim )) < MAX_ANGLE_ERROR

        I_residual_corr_mag_resid = (S1 - S2 + S3 - S4) / (4.0 * F__I_part_avg)
        I_residual_corr_hc_angle_1 = math.cos(I__I_part_avg_prelim)*delta_prelim*(alpha1-alpha2+(math.pi-alpha3)-(math.pi-alpha4))/4.0
        I_residual_corr_hc_angle_2 = math.sin( 2* I__I_part_avg_prelim)*( np.square(alpha1) + np.square(alpha2) + np.square(math.pi-alpha3) + np.square(math.pi-alpha4))/4.0

        I__I_part_avg = I__I_part_avg_ideal  - I_residual_corr_mag_resid - I_residual_corr_hc_angle_1 - I_residual_corr_hc_angle_2
        
        assert not np.isnan(I__I_part_avg)

        # if I__I_part_avg > math.pi/2.: I__I_part_avg = I__I_part_avg - math.pi

        I__I_part_avg_prelim = I__I_part_avg
        
        if mag_obs_dataset is not None:

            Z__I_part_avg__from_DIM_at_variometer = (F__I_part_avg - common_calc_params.F_pillar_diff_variometer) * np.sin(I__I_part_avg - common_calc_params.dI_pillar_diff) # DIM I rough variation transferred to variometer pillar
            # Other: Z__I_part_avg__from_DIM_at_variometer = np.mean( (F[ Incl_series_indeces ] - variometer_data_def.F_pillar_diff_variometer) * np.sin(I__I_part_avg - variometer_data_def.dI_pillar_diff)  )# DIM I rough variation trasnferred to variometer pillar

            Z_variometer_avg = np.mean(variometer_data['Z'][Incl_series_indeces])
            Z0 = Z__I_part_avg__from_DIM_at_variometer-Z_variometer_avg
            
            arcsin_parameter = (variometer_data['Z'] + Z0)/(F - common_calc_params.F_pillar_diff_variometer)
            
            try:  # np.arcsin(), if received values out of the valid range, normally gives a runtime warning only
              with np.errstate(invalid='raise'):
                I_var_at_variometer =  np.arcsin(arcsin_parameter)
            except FloatingPointError:
                # invalid_arcsin_parameters = np.abs(arcsin_parameter) > 1
                # print( "arcsin() error: ", np.flatnonzero(invalid_arcsin_parameters), arcsin_parameter[invalid_arcsin_parameters])
                # print("Time: ",  variometer_data[COMPNAME__TIMESTAMP][invalid_arcsin_parameters] ,"F: ", F[invalid_arcsin_parameters], "  Z: ", variometer_data['Z'][invalid_arcsin_parameters])
                # We get typically here if baseline corrected Z is higher than F, so arcsin receives invalid parameter
                return "Error calculating Inclination (baseline-corrected-Z > F)", None
            
            # 2026.05.04 experience: in a rare case, 'I_var_at_variometer' calculations may results in NANs!
            if np.any(np.isnan(I_var_at_variometer) ):
                return "I variation calculation error", None
            
            I_var_at_abs = I_var_at_variometer + common_calc_params.dI_pillar_diff
            
            epsilon_I_ideal  = np.pi/2 - np.mean((theta1,theta2,theta3,theta4)) # amdp.Epsilon_I(NU, SD, ND, SU)
            epsilon_I_residual_corr = (S1 - S2 - S3 + S4) / (4*F__I_part_avg)
            I1 = np.mean(I_var_at_abs [ NU_indeces ])
            I2 = np.mean(I_var_at_abs [ SD_indeces ])
            I3 = np.mean(I_var_at_abs [ ND_indeces ])
            I4 = np.mean(I_var_at_abs [ SU_indeces ])
            epsilon_I_var_corr = (I1+I2-I3-I4)/4.0
            epsilon_I  = epsilon_I_ideal + epsilon_I_residual_corr + epsilon_I_var_corr
    
            M_offset_I_ideal =  -((theta1-theta2-theta3+theta4)/4 + np.pi/2)*F__I_part_avg #   amdp.Offset_I(NU, SD, ND, SU, F__I_part_avg)
    
            # M_offset_I_var_corr =  (np.mean( I_var_at_abs[NU_indeces + ND_indeces] ) - np.mean( I_var_at_abs[SD_indeces + SU_indeces]))/2.0 *F__I_part_avg
            M_offset_I_var_corr =  (I1-I2+I3-I4)/4.0 *F__I_part_avg
            M_offset_I_residual_corr = (S1 + S2 + S3 + S4) / 4.0
            M_offset_I =  M_offset_I_ideal + M_offset_I_residual_corr + M_offset_I_var_corr
            
            I__D_part_avg  = np.mean( I_var_at_abs[Decl_series_indeces] )
            H__D_part_avg  = np.mean( F[Decl_series_indeces] * np.cos(I_var_at_abs[Decl_series_indeces]) ) 
            
            
        else:
            I__D_part_avg  = I__I_part_avg
            H__D_part_avg  = np.mean( F[Decl_series_indeces] * np.cos( I__D_part_avg  ) ) 

        # -----------------------------
        # Calculation of Declination
        # -----------------------------

        F__D_part_avg = np.mean( F[ Decl_series_indeces ] )

        S1 = np.mean( mag[EU_indeces] )
        S2 = np.mean( mag[WU_indeces] )
        S3 = np.mean( mag[ED_indeces] )
        S4 = np.mean( mag[WD_indeces] )

        theta1 =  angle_utils.select_equivalent_angle( vc_NED[EU_indeces], 0.0,   MAX_ANGLE_ERROR, False ) 
        if np.any(np.isnan(theta1)):
            # Maybe the reference angle is wrong??
            return "theta calculation error", None
        theta1 = np.mean( theta1 )
        theta2 = np.mean( angle_utils.select_equivalent_angle( vc_NED[WU_indeces], 0.0,   MAX_ANGLE_ERROR ) )
        theta3 = np.mean( angle_utils.select_equivalent_angle( vc_NED[ED_indeces], np.pi, MAX_ANGLE_ERROR ) )
        theta4 = np.mean( angle_utils.select_equivalent_angle( vc_NED[WD_indeces], np.pi, MAX_ANGLE_ERROR ) )


        # alpha1 = angle_utils.select_equivalent_angle( np.mean(hc_NED[EU_indeces]), np.pi/2.0,   MAX_ANGLE_ERROR )
        # alpha2 = angle_utils.select_equivalent_angle( np.mean(hc_NED[WU_indeces]), 3*np.pi/2.0, MAX_ANGLE_ERROR )
        # alpha3 = angle_utils.select_equivalent_angle( np.mean(hc_NED[ED_indeces]), 3*np.pi/2.0,   MAX_ANGLE_ERROR )
        # alpha4 = angle_utils.select_equivalent_angle( np.mean(hc_NED[WD_indeces]), np.pi/2.0,   MAX_ANGLE_ERROR )
        del alpha1,alpha2,alpha3,alpha4
        fi1 = np.mean(hc_NED[EU_indeces])
        fi2 = np.mean(hc_NED[WU_indeces])
        fi3 = np.mean(hc_NED[ED_indeces])
        fi4 = np.mean(hc_NED[WD_indeces])

        # New faster method tested from 2025.12-18
        H__D_part_avg_new = F__D_part_avg * np.cos(I__D_part_avg)
        if  abs(H__D_part_avg_new - H__D_part_avg)>0.02:
            return "H-check failed", None

        D__D_part_avg_ideal = np.mean((fi1,fi2,fi3,fi4)) - np.pi  # amdp.Declination(MIRA_WGS84, true_angle_calculator.last_measured_refpos_hc_theo_SU, true_angle_calculator.last_measured_refpos_hc_theo_SD, ED, WD, EU, WU)
        # Correction terms from Equation (12)
        D_residual_corr_mag = (S1 - S2 + S3 - S4) / ( 4.0 * H__D_part_avg)
        #if test_mode:
        #    D_residual_corr_mag = np.mean(S1/(F[EU_indeces] * np.cos(I_var_at_abs[EU_indeces]))  - S2/(F[WU_indeces] * np.cos(I_var_at_abs[WU_indeces])) + S3/(F[ED_indeces] * np.cos(I_var_at_abs[ED_indeces]))  - S4/(F[WD_indeces] * np.cos(I_var_at_abs[WD_indeces])))/4.0

        # D_residual_corr_vc_angle =  np.tan(I__D_part_avg )*( angle_utils.angle_distance(theta2,theta1) + angle_utils.angle_distance(theta3,theta4) )/ 4.0
        D_residual_corr_vc_angle =  np.tan(I__D_part_avg) *( theta1-theta2-theta3+theta4 )/ 4.0 # Note: in eq.12. of Lauridsen's paper equation, the "theta" part is included in the tangent function.
                    # However, test runs with simulated dat indicates, that this part must be outside from the tangent function.

        # if test_mode:
        #    D_residual_corr_vc_angle = 0
        D__D_part_avg = D__D_part_avg_ideal + D_residual_corr_mag - D_residual_corr_vc_angle

        if D__D_part_avg > math.pi: D__D_part_avg =  D__D_part_avg - 2 * math.pi

        if mag_obs_dataset is not None:
            # Calculation of baselines
            
            DIM_HoizTot__I_port_avg__at_variometer = (F__I_part_avg - common_calc_params.F_pillar_diff_variometer) * np.cos(I__I_part_avg - common_calc_params.dI_pillar_diff)        
            if baseline_calc_approach == "D0-angle":
                H0 = np.sqrt( np.square(DIM_HoizTot__I_port_avg__at_variometer) - np.square( np.mean(variometer_data[ D_comp_name  ][Incl_series_indeces])  )  )  - np.mean(variometer_data[ H_comp_name ][Incl_series_indeces])
                D_variometer_without_D_baseline = np.mean(np.arctan2( variometer_data[D_comp_name][Decl_series_indeces],  variometer_data[H_comp_name][Decl_series_indeces] + H0 ))
                D0 = D__D_part_avg - common_calc_params.dD_pillar_diff - D_variometer_without_D_baseline
                D_var_at_variometer =  np.arctan2( variometer_data[D_comp_name],  variometer_data[H_comp_name] + H0) + D0            
            elif baseline_calc_approach == "Y0-nT":
                DIM_X__at_variometer = DIM_HoizTot__I_port_avg__at_variometer * np.cos(D__I_part_avg_prelim - mag_obs_dataset.get_decbas())
                DIM_HorizTot__D_part_avg__at_variometer = (F__D_part_avg - common_calc_params.F_pillar_diff_variometer) * np.cos(I__D_part_avg - common_calc_params.dI_pillar_diff)                    
                DIM_Y__at_variometer = DIM_HorizTot__D_part_avg__at_variometer * np.sin(D__D_part_avg - mag_obs_dataset.get_decbas() )
                H0 = DIM_X__at_variometer - np.mean(variometer_data[ H_comp_name ][Incl_series_indeces])
                #if iteration_count==MAX_ITERATIONS-1:
                #    pass
                D0  = DIM_Y__at_variometer  - np.mean(variometer_data[ D_comp_name ][Decl_series_indeces])
                # ToDo: D0 számítható abból is, hogy HorizTotal = X^2  + Y^2 
                D_var_at_variometer = np.arctan2( variometer_data[D_comp_name] + D0,  variometer_data[H_comp_name] + H0) + mag_obs_dataset.get_decbas()
            else:
                assert False
            
            D_var_at_abs = D_var_at_variometer + common_calc_params.dD_pillar_diff

            # Equation (13) 'D1+D2-D3-D4' part:
    
            D1 = np.mean( D_var_at_abs[EU_indeces] )
            D2 = np.mean( D_var_at_abs[WU_indeces] )
            D3 = np.mean( D_var_at_abs[ED_indeces] )
            D4 = np.mean( D_var_at_abs[WD_indeces] )
    
            # Equation (13) delta
            # delta_D_ideal  = amdp.Delta(ED, WD, EU, WU)  # QQQ: nem fordított az előjel????
            # delta_D_corr_1 = (D1+D2-D3-D4)/4
            # delta_D_corr_1 = -delta_D_corr_1 # QQQ ellenkezo elojel kell!
            # delta_D_corr_2 = (S1-S2-S3+S4)/H__D_part_avg/4
            # delta_D_corr_3 = angle_utils.angle_distance(theta2,theta1) + angle_utils.angle_distance(theta4,theta3) /4.0 * math.tan(I__D_part_avg)
            # delta_D   = delta_D_ideal +    delta_D_corr_1 + delta_D_corr_2 + delta_D_corr_3
            delta_D = 0.25 * (D1+D2-D3-D4-fi1-fi2+fi3+fi4 - (S1-S2-S3+S4)/H__D_part_avg) + (theta1-theta2+theta3-theta4)/4.0 * math.tan(I__D_part_avg)
    
            # epsilon_D_ideal  = amdp.Epsilon_D(ED, WD, EU, WU, I__D_part_avg)
            # epsilon_D_corr_1 = (-D1+D2+D3-D4)         / ( 4* math.tan(I__D_part_avg) )  # !!! tapasztalat szerint itt minusz kell!!!
            # H__D_part_avg = math.cos(I__D_part_avg)  * F__D_part_avg
            # epsilon_D_corr_2 = (S1+S2-S3-S4)/H__D_part_avg/( 4* math.tan(I__D_part_avg) )
            # epsilon_D = epsilon_D_ideal + epsilon_D_corr_1  + epsilon_D_corr_2
    
    
            # epsilon_D is affected by hc eccentricity error!!!!
            epsilon_D = 1.0 / (4.0*math.tan(I__D_part_avg)) * ( fi1 - fi2 - fi3 + fi4 - D1 + D2 + D3 - D4 + (S1+S2-S3-S4)/H__D_part_avg + 2 * np.pi)  + np.pi/2 - np.mean((theta1,theta2,theta3,theta4))
            #if True:
            #    epsilon_D_1 = (- D1 + D2 + D3 - D4) / (4.0*math.tan(I__D_part_avg))
            #    epsilon_D_2 = (S1+S2-S3-S4)/H__D_part_avg / (4.0*math.tan(I__D_part_avg))
            #    epsilon_D_3 = np.mean((theta1,theta2,theta3,theta4))
    
            #    print( ">>> {:.3f} {:.3f} {:.3f}".format(math.degrees(epsilon_D_1), math.degrees(epsilon_D_2), 90-math.degrees(epsilon_D_3) ))
    
            # M_offset_D = amdp.Offset_D(ED, WD, EU, WU, H__D_part_avg)  # Negative sign added at 2024.03.11, then removed at 2024.03.18
            # # H_avg_at_decl_series = F__D_part_avg * np.cos(I__I_part_avg)
            # M_offset_D = M_offset_D - ((np.mean( D_var_at_abs[ED_indeces+EU_indeces]) - np.mean( D_var_at_abs[WD_indeces + WU_indeces]))/2.0 * H__D_part_avg)
            # M_offset_D = M_offset_D + (S1+S2+S3+S4)/4 - ((theta1+theta2-theta3-theta4)/4 + math.pi/2) * (Z_variometer_avg + Z0) # QQQ: Z értéke!!!
    
            Z_variometer_avg_D = np.mean(variometer_data['Z'][Decl_series_indeces]) + Z0
            M_offset_D = (fi1-fi2+fi3-fi4)/4.0 * H__D_part_avg - (D1-D2+D3-D4)/4.0 * H__D_part_avg + (S1+S2+S3+S4)/4.0 - (( theta1+theta2-theta3-theta4 ) / 4.0 + np.pi/2.0) * (Z_variometer_avg_D) # QQQ: Z értéke!!!
    
            # delta_D     -= (np.mean( D_var_at_abs[ED_indeces + WD_indeces]) - np.mean(D_var_at_abs[EU_indeces + WU_indeces]))/2.0 # QQQ
            # delta     -=
            # epsilon    = np.mean( [ epsilon_D, epsilon_I]  )
    
            # However epsilon_D is affected by hc eccentricity error, in the full equation system, it seems to be better to use the average!!!!
            delta        = delta_D
            delta_prelim = delta_D
    
            D__I_part_avg = np.mean( D_var_at_abs[Incl_series_indeces] )
            D__I_part_avg_prelim = D__I_part_avg
    
            # 'Decl_series_indeces' and 'Incl_series_indeces' holds all raw points.
    
    
            # Determining the scale.
            # Note: offset depends on scale
            # mean_offset = np.mean([M_offset_I,M_offset_D])

        if calc_params.separate_D_I_processing:
            # Szimuláicók szerint ez ugyanolyan végeredményt ad, de a meas points resid ennél a számításnál simább, ha van hc eccentr hiba!
            M_offset_common = np.empty( (len(F),) )
            epsilon_common  = np.empty( (len(F),) )
            M_offset_common[Decl_series_indeces] = M_offset_D
            M_offset_common[Incl_series_indeces + SC_indeces] = M_offset_I
            epsilon_common[Decl_series_indeces] = epsilon_D
            epsilon_common[Incl_series_indeces + SC_indeces] = epsilon_I
        else:
            M_offset_common    = (M_offset_I+M_offset_D)/2.0
            epsilon_common    =  (epsilon_D + epsilon_I)/2.0

        if b_estimate_M_scale: M_offset_common = M_offset_common * M_new_scale # Tapasztalat szerint ez gyakrolatilag nem változtat semmit

        M_with_unit_scale = amdp.dim_instrument_model(F, vc_NED,  hc_NED, I_var_at_abs, D_var_at_abs, epsilon_common, delta, M_offset_common, s = 1.0)
        
        if np.any(np.isnan(M_with_unit_scale)):
            assert False
            
        if b_estimate_M_scale:
            remaining_offset, M_new_scale, resid_info = mu.LSQ( M_with_unit_scale, measured_mag  )
            
            if debug_iteration_process:
                iteration_progress_data[iteration_count]['I_I'] = math.degrees(I__I_part_avg)
                iteration_progress_data[iteration_count]['D_I'] = math.degrees(D__I_part_avg)
                iteration_progress_data[iteration_count]['delta_arcsec'] = angle_utils.rad_to_arcsec(delta_D)
                iteration_progress_data[iteration_count]['M_scale'] = M_new_scale
    
            # M_new_scale = np.average(  mag[indeces_used_for_scale_estimation] / M_with_unit_scale[indeces_used_for_scale_estimation], weights = np.abs(measured_mag[indeces_used_for_scale_estimation]) )
    
            #if (abs(M_new_scale)<0.2) or (abs(M_new_scale)>1.8):
            #    return "M_scale calculation error", None
    
            # print(1*(M_new_scale if test_mode else 1.0), remaining_offset, M_new_scale,iteration_count)
            # TODO if (abs(M_new_scale)-1.0 < 1e-3) and (iteration_count>2) : break
    
            # sc += [M_new_scale]
            # if test_mode and  (iteration_count==1): break
            if calc_params.optimize_M_scale:
                mag = measured_mag / M_new_scale

    M_with_unit_scale = amdp.dim_instrument_model(F, vc_NED,  hc_NED, I_var_at_abs, D_var_at_abs, epsilon_common, delta, M_offset_common, s = 1.0)
    if b_estimate_M_scale:
        remaining_offset,M_new_scale, resid_info = mu.LSQ( M_with_unit_scale, measured_mag  )
        res_mag_sensor.scale =  M_new_scale # np.mean(measured_mag[indeces_used_for_scale_estimation] / mag[indeces_used_for_scale_estimation])
    else:
        res_mag_sensor.scale =  M_new_scale

    if debug_iteration_process:
        plot_utils.plot_xy( range(len(MAX_ITERATIONS)), iteration_progress_data )
    # input("Press Enter to continue...")

    # if abs( abs(res_mag_sensor.scale)-1.0) > 0.5:
    #    assert False, f"{res_mag_sensor.scale=}"
    res_mag_sensor.delta_arcsec   = angle_utils.rad_to_arcsec(delta)
    res_mag_sensor.epsilon_arcmin = angle_utils.rad_to_arcmin(epsilon_I) # epsilon_D is affected by hc eccentricity error

    M_offset_I = M_offset_I *  res_mag_sensor.scale
    M_offset_D = M_offset_D *  res_mag_sensor.scale
    res_mag_sensor.offset  = np.mean([M_offset_I,M_offset_D])

    # Calculating residuals:
    calc_mag = amdp.dim_instrument_model(F, vc_NED, hc_NED, I_var_at_abs, D_var_at_abs, epsilon_common, delta, M_offset_common, res_mag_sensor.scale)
    resid = measured_mag - calc_mag
    resid_std = np.std(resid)

    decl_incl_part_separate_results = \
    [
       { "name":"M-offset", "unit":"nT",     "value-D":  M_offset_D,                   "value-I": M_offset_I,                   "range_min":res_mag_sensor.offset-3, "range_max":res_mag_sensor.offset+3},
       { "name":"epsilon",  "unit":"arcmin", "value-D":  angle_utils.rad_to_arcmin(epsilon_D),  "value-I": angle_utils.rad_to_arcmin(epsilon_I),  "range_min":None,       "range_max":None},
    ]
    
    return None, (I_var_at_abs, D_var_at_abs, I__I_part_avg, I__D_part_avg, D__I_part_avg, D__D_part_avg, Incl_series_indeces, Decl_series_indeces, F__I_part_avg, F__D_part_avg, resid_std, resid, decl_incl_part_separate_results, H0, D0, Z0)



# This function calculates the pseudo 0 nT positions of a residual measurement.
# Recently the calculation is based on linear interpolation of more measurements
# in each main position: this measureing method is used by ADS.
# hc/vc angles can be radans and degrees as well.
# Offset/Delta/epsilon is not taken into account, so 0 nT positions are the positions
# where the DIM would measures 0nT, and not where the field is really zero.
# Similar function in deprecated code: 'NM_calculate_zeropos_angles_and_positions'
def calc_null_method_positions(NM_position_indeces, M, hc_theo, vc_theo, data_for_time_reduction):
                
    rotated_angle_null_positions = {}
    rotated_angle_null_positions_zeropos_residstd = {}
    fixed_angle_null_positions = {}
    
    # We keep the original measurement order in the returned dict, so we expect time ordered list!
    
    for pos in NM_position_indeces:
        
        pos_indeces = NM_position_indeces[pos]
        if len(pos_indeces) < 2: continue # Scale check may have one readings        
        if pos[0] in ['E','W']:  # East/West positions of declination series
            # Declination series
            fixed_angles   = vc_theo
            rotated_angles = hc_theo
        elif pos[0] in ['N','S']:  # North/South positions of inclination series
            fixed_angles   = hc_theo
            rotated_angles = vc_theo
        else:
            assert False

        residstd_list = ""
        for corr in [0]:  #  [ -1, 0, +1 ]:

            if data_for_time_reduction is not None:
                hc_NED = data_for_time_reduction["hc_NED"][pos_indeces]
                vc_NED = data_for_time_reduction["vc_NED"][pos_indeces]
                abs_magfield = gu.magfield_class( I = data_for_time_reduction["magfield"].I[pos_indeces],  D = data_for_time_reduction["magfield"].D[pos_indeces], F = data_for_time_reduction["magfield"].F[pos_indeces] )
                sensor = data_for_time_reduction["mag-sensor"]
                M_time_reducted = gu.M_recordings_time_reduction( 0 * len(pos_indeces)//2 ,  M[pos_indeces], hc_NED, vc_NED, abs_magfield, sensor, correction = corr )
            else:
                M_time_reducted = M[pos_indeces]
            
            fixed_angle_null_positions[pos] = np.mean( fixed_angles[pos_indeces] ) 
            offset,scale,[residstd, resids] = mu.LSQ( rotated_angles[pos_indeces], M_time_reducted, return_resid_stddev = True, return_resid_array = True  )
            null_pos = -offset/scale
            rotated_angle_null_positions[pos] = null_pos        
            residstd_list += "{:}:{:.5f} ".format(corr, residstd)
        
        rotated_angle_null_positions_zeropos_residstd[pos] = residstd_list
        
        
        if False and ( (null_pos < np.min(rotated_angles[pos_indeces])) or (null_pos > np.max(rotated_angles[pos_indeces])) ):
            # print( f"{null_pos}  {rotated_angles[pos_indeces]}")
            plot_utils.plot_xy(rotated_angles[pos_indeces], M_time_reducted, "rotated-angle", "M [nT]","o", close_plot = False) 
            fitted_line_x = np.linspace( min(rotated_angles[pos_indeces]), max(rotated_angles[pos_indeces]), 10, endpoint = True)
            fitted_line_y = offset + scale * fitted_line_x
            plot_utils.plot_xy(fitted_line_x, fitted_line_y, "rotated-angle", "M [nT]","-", close_plot = False) 
            fitted_line_x = np.array([null_pos])
            fitted_line_y = offset + scale * fitted_line_x
            plot_utils.plot_xy(fitted_line_x, fitted_line_y, "rotated-angle", "M [nT]","x", close_plot = True) 
            
            plot_utils.plot_xy(range( len(resids) ), resids, "rotated-angle", "M-resid [nT]","x", close_plot = True)             
    
            
        
    # "SC" stands for scale-check
    SC_indeces = NM_position_indeces[du.SCALECHECK_POS_ID]
    if len(SC_indeces) == 0: 
        # There was no separate scale check measurement => Selecting the highest vc angles's from the South/Face-LeftL ("SU") measurement
        SC_pos =  "SL"
        pos_indeces = NM_position_indeces[SC_pos]
        assert len(pos_indeces) > 1
        highest_vc_index = np.argmax(vc_theo[pos_indeces])
        SC_vc = vc_theo[pos_indeces][highest_vc_index]
        SC_hc = hc_theo[pos_indeces][highest_vc_index]
        SC_M  = M[pos_indeces][highest_vc_index]
    else:
        assert len(SC_indeces) == 1
        # There was one scale measurement
        SC_vc = vc_theo[SC_indeces[0]]
        SC_hc = hc_theo[SC_indeces[0]]
        SC_M  = M[SC_indeces[0]]
        
    return rotated_angle_null_positions, fixed_angle_null_positions, [SC_vc, SC_hc, SC_M], rotated_angle_null_positions_zeropos_residstd



# --------------------------------------------------------------------------------------
# This main function runs a separate optimization for each absolute measurements (referred here as processing-windows), 
# to get spot baseline values and sensor offset/scale and misaligment results.
# All parameters (M_scale/M_offset/epsilon/delta...) are re-calcualted for each processing window.
# This is the conventional data processing/evaluation method for ADS absolute measurements.
# Data processing options are defined in the 'calc_params_list' array: it can contain one or more option set.
# Optional processing type is the "all-at-once" optimization, that is temporarly not available.
# --------------------------------------------------------------------------------------
def dim_data_process(dat_files_format, dat_files_path, dim_datafile_name_def, meas_pos, DIM_instr, mag_obs_dataset, calc_params_common, calc_params_list, ui_preferences, processing_mode = "PER-CYCLES", aux_data = None):

    assert isinstance(DIM_instr,dim_instr.DIM_instrument_model_class)    
    assert (mag_obs_dataset is None)     or  isinstance(mag_obs_dataset, geomag_obs_dataset.mag_obs_dataset)    
    assert isinstance(calc_params_list,list)
    assert isinstance(calc_params_list[0], spot_baselines_calculation_parameters_class)
    assert processing_mode in ["PER-CYCLES", "ALL-IN-ONE"]
    assert isinstance(meas_pos, gu.measurement_pos_class ) # or (meas_pos is None) 

    # -------------------
    # -- UI parameters --
    # -------------------
    # Some variables to enable/disable some user-interface behaviour:
    input_data_plots    = ui_preferences.pop("input-data-plots", False)
    enable_cycle_plots  = ui_preferences.pop("processing-window-plots", False)
    enable_cycle_tables = ui_preferences.pop("processing-window-tables", False)
    enable_post_plots   = ui_preferences.pop("post-plots", False)
    cu.wait_for_enter_at_infos = ui_preferences.pop("wait-for-enter-at-infos", False)
    cu.wait_for_enter_at_warnings = ui_preferences.pop("wait-for-enter-at-warnings", False)
    
    immediately_plot_1D_residual_for_each_cycle = False # enable_cycle_plots  => it is plotted after table is printed out
    # do_corr_variation_plot             = enable_cycle_plots and (varcorr_data_def is not None) # Active only if there is ref -data
    do_2D_residual_plot_for_each_cycle = False and enable_cycle_plots # residual in hc/vc plane
    
    
    # ---------------------------------
    # -- Rarerly changed parameters --
    # ---------------------------------
    glob_test_mode = ""
    do_print_new_random_step_results = False
    do_tick_mod_4_analysis = (DIM_instr.type == "ADS") and ("do_tick_mod_4_analysis" in calc_params_common.post_processing_steps)

    # --------------------------
    # --- End of parameters  ---
    # --------------------------
    
    warn_logs = []
    saved_files = []
    
    DIM_meas_cycle_magfield__common = gu.magfield_class(meas_pos)
    nof_calc_params = len(calc_params_list)    
    b_more_calc_params = nof_calc_params > 1

    # Checking some settings of calc. variations and determining which magnetoemter sensors we need to load.
    mag_indeces = []
    mag_sensor_names = []
    for calc_params in calc_params_list:
        if not (calc_params.mag_sensor_index in mag_indeces):
            mag_indeces = mag_indeces + [calc_params.mag_sensor_index]
            assert calc_params.mag_sensor_index == 0  # Support for 3 component sensor is temporarly unavailable!
            
            mag_sensor_names = [ "M0" ]
            
            # if calc_params.mag_sensor_index > max_mag_index :  max_mag_index = calc_params.mag_sensor_index
            # if calc_params.mag_sensor_index==0:
            #    mag_sensor_names = mag_sensor_names + [ dim_instr.def_mag_sensor_name ]
            # else:
            #    mag_sensor_names = mag_sensor_names + ["MAG-2" ]

    DIM_data_loader = dim_measurement_loader.DIM_dat_file_reader_class()
    
    dim_datafile_name_def = gu.replace_masks_in_filename( dim_datafile_name_def, meas_pos )
        
    
    DIM_data_loader.load(dat_files_format, dat_files_path, dim_datafile_name_def, DIM_instr, calc_params_common, calc_params_list, meas_pos)
    # Note: 'DIM_data_loader' will fill the vc0 field of 'DIM_instr' in case of ADS raw data file readins

    cu.wait_for_enter()
    print("")  # Adding a line after DIM data files' plots

    print("DIM instrument's parameters:\n\tname: {:}\n\tmag-sensor(s): {:}\n\tvc0: {:.3f}°".format(DIM_instr.name, mag_sensor_names, math.degrees(DIM_instr.vc0)))
    print()
    dim_instr.last_applied_vc0 = DIM_instr.vc0  # It is ' DIM_data_loader.vc0_from_file '

    if calc_params_common.do_self_checks:
        ref_vc_angle = DIM_instr.get_horizontal_vc_angle(DIM_instr.hc_angle_reference_sensor_pos(),"rad")
        avc = (dim_azimuth.v_theo_to_ned( ref_vc_angle  ) + math.tau) % math.tau
        assert ( math.degrees(avc)<1.0)

    cu.wait_for_enter()

    # Preparing processing windows (for each calc. var):
    nof_concatenated_file_cycles = 1 # dim_measurement_loader.ALL_MEASUREMENTS_TOGETHER if  processing_mode == "ALL-IN-ONE" else 1
    DIM_data_loader.prepare_processing_windows(nof_concatenated_file_cycles)
    cu.wait_for_enter()

    # Plotting input file data (with calculation processing-window limits)
    
    if input_data_plots: DIM_data_loader.plot_load_data()
    DIM_data_loader.plot_mag_meas_slopes(None,None) # QQQ Plotting a preliminary, non-variometer corrected slopes
    cu.wait_for_enter()
    

    def_magfield = gu.magfield_class( meas_pos ) # 'meas_pos' is added in 2025.11.14
    DIM_data_loader.check_load_data( def_magfield )
    cu.wait_for_enter()

    # Loading variometer data
    if (mag_obs_dataset is None) or (not mag_obs_dataset.has_variometer_data() ):
        variometer_data_parallel_mag_data_index = -1
    else:
        DIM_data_loader.load_parallel_mag_recordings(mag_obs_dataset, time_shift_sec = calc_params_common.variometer_time_shift_sec, use_subsecond_interpol = (calc_params_common.DIM_time_shift_sec is not None) )
        cu.wait_for_enter()
        
    has_valid_window = DIM_data_loader.print_processing_windows()
    if not has_valid_window:
        warn_logs += DIM_data_loader.get_warnings()
        warn_logs += [ "There is no calculation window that could be processed - no result files saved!" ]
        return None, None, None, warn_logs

    cu.wait_for_enter()
        

    #if (mag_obs_dataset is None) or (not mag_obs_dataset.has_F_data() ) or (not isinstance(mag_obs_dataset.F_data, geomag_obs_dataset.mag_obs_data_reader)):
    #    F_data_parallel_mag_data_index = -1
    #else:
    #    F_data_parallel_mag_data_index = DIM_data_loader.load_parallel_mag_recordings(mag_obs_dataset.F_data, "F-data", time_shift_sec = F_time_shift_sec)


    print("")  # Adding a line after ADS data plots

    #  'mag_obs_dataset' knows its frame after data is load
    for calc_param in calc_params_list:
        calc_param.prepare_for_calc(mag_obs_dataset)
        
    print()
    print("------ Calculating spot baseline points ------")


       
    
    # Record for point-wise results in 'ra_measpoints' list of arrays.
    # Each calc.var. has its own ndarray inside.
    # At start, we reserving memory for all input raw records: later this array may be shrinked.
    angle_unit = "deg"
    ra_measpoints_comps = [
        ("datetime",        "f8", f"%.{mu.TIMELABEL_DECIMAL_DIGIT_FOR_SECOND}f" ),
        ("win_nr",          "u4", "%u"), 
        ("win_sample_nr",   "u4", "%u"), 
        ("pos_ID",          "u4", "%u"),
        ("hc_THEO",         "f8", f"%.{angle_utils.angle_expression_decimal_digits[angle_unit]}f", angle_unit  ),
        ("vc_THEO",         "f8", f"%.{angle_utils.angle_expression_decimal_digits[angle_unit]}f", angle_unit  ),
        ("hc_NED",          "f8", f"%.{angle_utils.angle_expression_decimal_digits[angle_unit]}f", angle_unit  ),
        ("vc_NED",          "f8", f"%.{angle_utils.angle_expression_decimal_digits[angle_unit]}f", angle_unit  ),
        ("M",               "f8", "%.2f",  "nT"), 
        ("M_raw_resid_nT",  "f8", "%.3f",  "nT"), 
        ("M_mean_removed_resid_nT", "f8",  "%.3f", "nT"),
        # ("epsilon","f8", 2, arcmin), ("delta","f8", 1, "arcsec"), ("M_offset","f8", 2, nT), ("M_scale","f8", 3, "nT/nT"),
        ("M_meas_noise",'f8', "%.2f", "nT"), ("M_meas_slope",'f8', "%.2f", "nT/nT"), 
        # ("M_slope_from_variometer",'f8'),
        # ("varcorr_nT",float),
        #("F",        "f8", 2,  "nT"),
        #("I","f8"), ("D","f8"), 
        # ("var_HX",float), ("var_DY",float), ("var_Z",float) 
        ]
    
    ra_measpoints_dtypes = [ item[:2] for item in ra_measpoints_comps ]
    ra_measpoints_header = [] 
    ra_measpoints_fmt = []
    for comp_def in ra_measpoints_comps:
        ra_measpoints_fmt += [ comp_def[2] ]
        if len(comp_def) < 4 or (comp_def[3] == ""):
            ra_measpoints_header += [ comp_def[0] ]
        else:
            ra_measpoints_header += [ comp_def[0]+"["+comp_def[3]+"]" ]
    ra_measpoints_header = " ".join(ra_measpoints_header)  # Converting list to string
    
    
    if do_tick_mod_4_analysis: 
        ra_measpoints_dtypes += [ ("hc_tick_mod_4", 'u4', 0, ""), ("vc_tick_mod_4", 'u4', 0, "") ]
   

    XH_compname = mag_obs_dataset.usage_frame[0]
    YD_compname = mag_obs_dataset.usage_frame[1]
    
    result_array_dtypes, result_array_comps_params, concatenated_plot_y_components = get_amdp_results_array_data(calc_params_common.baseline_calc_approach, XH_compname, YD_compname, meas_pos.pillars, DIM_data_loader.get_aux_comp_names(True), calc_params)   
    for icalcvar, calcvar in enumerate(calc_params_list):
        calcvar.result_array = np.empty( DIM_data_loader.number_of_processing_windows, dtype = result_array_dtypes) #  mu.datarecord_array_class(DIM_data_loader.number_of_processing_windows) 
        # Emptying result array:
        for comp_name in calcvar.result_array.dtype.names:
            calcvar.result_array[comp_name][:] = result_array_comps_params[comp_name].empty_value
        calcvar.ra_measpoints = np.zeros( ( DIM_data_loader.nr_of_load_records() ), dtype = ra_measpoints_dtypes    ) 
        for comp_name in calcvar.ra_measpoints.dtype.names:
            if calcvar.ra_measpoints.dtype[comp_name].kind == 'f':   # 'f' = floating point
                calcvar.ra_measpoints[comp_name] = np.nan        
        calcvar.nr_of_ra_measpoints_records = 0
        calcvar.varcorr_prev_M_scale =  None 
        
        if calcvar.calculated_null_positions_file is not None:
            calcvar.calculated_null_positions = {}
        
        
        
    # remove_mean_comps = ['H0', 'X0', 'D0', 'D0-angle','Y0', 'D0_nT_from_angle','Z0', 'M_offset', 'M_scale', 'epsilon', 'delta' ] if calc_params_common.b_remove_means else []
    remove_mean_comps = [ "azimark-hc-change_mean-removed"  ]
    remove_mean_comps += ["M-offset"]
    
    if calc_params_common.remove_means is not None:
        for remove_keyword in calc_params_common.remove_means:
            if type(remove_keyword) is bool:
                assert not remove_keyword  # "remove_means = no" results in a "False" setting, that is then put in a list of "[ False ]"
                assert len(calc_params_common.remove_means) == 1
            else:
                match remove_keyword:
                    case "baselines":
                            remove_mean_comps += baseline_components
                    case "no":
                        break
                    case _:  
                        assert False, f"{remove_keyword=}"

    continuous_mode = "waitforenter"

    # resid_per_pos_id_avgsum = [ ]   # Calculating average residual value for each pos-id
    # resid_per_pos_id_count  = [ ]   # Calculating average residual value for each pos-id
    for calc_params in calc_params_list:
        calc_params.resid_per_pos_id_avgsum = dict()  # Calculating average residual value for each pos-id
        calc_params.resid_per_pos_id_count  = dict()  # Calculating average residual value for each pos-id
 
    for processing_window_index in range(DIM_data_loader.number_of_processing_windows):
        
        # Used until 2026.09.16: common variaton plot for all calc. vars
        # if enable_cycle_plots:
        #     mu.print_text_in_frame("Processing window {:}{}".format(processing_window_index,  "- Data for all-calcvars' period"  if b_more_calc_params else "" ),None,"===")
        #     # Collecting togheter all calc vars's mag  record indeces, to have all timelabels
        #     mag_meas_file_records_indeces_all_calcvars = []
        #     all_calvcars_first_timelabel, all_calvcars_last_timelabel = None, None
        #     for calc_params_index, calc_params in enumerate( calc_params_list ):
        #         window_data = calc_params.processing_window_data[processing_window_index]
        #         mag_file_record_indeces = window_data.M_record_indeces
        #         mag_meas_file_records_indeces_all_calcvars = list( set(mag_meas_file_records_indeces_all_calcvars) | set(mag_file_record_indeces)   ) 
        #         if calc_params_index == 0:
        #             all_calvcars_first_timelabel, all_calvcars_last_timelabel = window_data.first_M_timelabel,     window_data.last_M_timelabel
        #         else:
        #             all_calvcars_first_timelabel = min (all_calvcars_first_timelabel, window_data.first_M_timelabel)
        #             all_calvcars_last_timelabel  = max (all_calvcars_last_timelabel, window_data.last_M_timelabel)
                    
        #     mag_meas_point_datetimes = DIM_data_loader.M_readings[ mag_meas_file_records_indeces_all_calcvars ] [ dim_measurement_loader.COMPNAME__TIMESTAMP ]
        #     mag_obs_dataset.plot_data_for_date_range( all_calvcars_first_timelabel, all_calvcars_last_timelabel, mark_datetimes = mag_meas_point_datetimes)

        for calc_params_index, calc_params in enumerate( calc_params_list ) :

            processing_window_data, meas_method_info, timelabel, mags, F_abs_pillar, hc_raw_tick, vc_raw_tick, hc_theo_deg, vc_theo_deg, hc_theo_rad, vc_theo_rad, hc_NED, vc_NED, meas_time, temperatures, processing_window_raw_data = DIM_data_loader.processing_cycle_data_select(DIM_instr, processing_window_index, calc_params_index, mag_indeces = mag_indeces, mag_precorr_functions = calc_params.mag_precorr_funcs)
            # cycle_name   = processing_window_data.name
            cycle_pillar = processing_window_data.pillar_name
            cycle_file_id_info = processing_window_data.file_id_str

            mid_datetime = (timelabel[0] + timelabel[-1]) / 2
            n = len(timelabel)

            if aux_data is not None:
                time_diff = abs(aux_data["datetime"] - mid_datetime)
                aux_data_idx_for_this_cycle = time_diff.argmin()
                assert (aux_data_idx_for_this_cycle>0) and (aux_data_idx_for_this_cycle < time_diff.shape[0]), f"{aux_data_idx_for_this_cycle=}"

            datetime_Ipart = datetime_Dpart = mid_datetime  # Can be overwritten later
            
            
            def add_result(comp_name, value, starting_val = None, range_min = None, range_max = None, ):
                nonlocal calc_params
                if isinstance(value, (np.bytes_, bytearray) ): value = value.decode()
                if result_array_comps_params[comp_name].add_to_table:
                    if isinstance(value,float):
                        mask="{:."+str(result_array_comps_params[comp_name].decimal_digits)+"f}"
                    else:
                        mask = "{}"
                    if result_array_comps_params[comp_name].optimized_param: 
                        comp_name_in_table = comp_name + " *"
                    else:
                        comp_name_in_table = comp_name
                    assert (range_min is None)  == (range_max is None)
                    table.add(comp_name_in_table, result_array_comps_params[comp_name].unit, mask, value, 1.0, starting_val, range_min, range_max)
                    
                # if result_array_comps_params[comp_name].add_to_ra:
                calc_params.result_array[processing_window_index][comp_name] = value 
            
            
            add_result("error_code",        CALC_ERROR_CODE__MISC_ERROR)
            add_result("program_type_idx",  meas_method_info[0]) 
            add_result("win_nr",            processing_window_index)  
            # add_result("win_name",          cycle_name )
            add_result("pier",              cycle_pillar )
            add_result("input_file_win_idx",      cycle_file_id_info ) 
            timestamp = mu.exceltime_2_datetimestr(timelabel[0])
            timestamp = timestamp.replace(' ', '_').replace('\t', '_')
            add_result("datetime-start-str",timestamp)
            add_result("datetime-start",    timelabel[0])  
            add_result("datetime-end",      timelabel[-1]) 
            add_result("datetime-mid",      mid_datetime) 
            add_result("datetime-Incl",     datetime_Ipart)  # Can be overwritten later
            add_result("datetime-Decl",     datetime_Dpart)  # Can be overwritten later
            # hour_of_day = (mid_datetime - int(mid_datetime))*24.0
            # add_result("hours-of-day",    hour_of_day)
            

            for temp_ch_name in temperatures:  # 'temperatures' is a dictionary
                # "T_electronics" / "T_room"
                temp_data = temperatures[temp_ch_name]
                add_result(temp_ch_name,   np.average(temp_data) )
                # add_result_comp(False,True, temp_ch_name+"-change",  "°C",  np.average(temp_data[-2:]) - np.average(temp_data[:2]), comp_prec = 1, plot = True  )
                
            # Getting out the values of the selected magnetometer sensor (needs for multi-magnetometer sensor configuration tests)
            dim_mag = mags[calc_params.mag_sensor_index].copy()
            
            calc_params_info = "" if not b_more_calc_params else "CALC.VAR:\"{}\"".format(calc_params.name)
            DIM_data_loader.print_processing_windows_data_info(processing_window_data,  enable_cycle_plots or enable_cycle_tables, first_line_extra_str="", calc_params_info = calc_params_info )
            
            cycle_skip_reason = None
            
            def cycle_check(cond, errcode, errmsg):
                nonlocal cycle_skip_reason
                if cond: return True
                cycle_skip_reason = errmsg
                calc_params.result_array[processing_window_index]["error_code"] = errcode 
                return False
            
            for do_calculation in [1]: # Dummy cycle, from which we can exit from a 'break' statement any time            
            
                for cycle_filter_item in calc_params.cycle_filters:
                    if cycle_filter_item[0] == CYCLE_FILTER__PILLAR_NAME:
                        if not cycle_check(cycle_pillar == cycle_filter_item[1], CALC_ERROR_CODE__PILLAR_MISMATCH, "Cycle is measured on other pillar {:} than {:}".format(cycle_pillar,cycle_filter_item[1])): break
                if cycle_skip_reason is not None: break
            
                if not cycle_check( not processing_window_data.excluded_by_data_selection_filter, CALC_ERROR_CODE__EXCLUDED_BY_FILTER,  "Window is excluded by a data selection filter"):
                    break
            
                if not cycle_check( processing_window_data.complete_azi_meas_data is not None, CALC_ERROR_CODE__NO_VALID_REFMARK, "Window has no valid reference mark reading data"): 
                    break

                if not cycle_check( not processing_window_data.missing_variometer_data, CALC_ERROR_CODE__MISSING_VARDATA, "Missing variometer data"): 
                    break
                                   
                if not cycle_check( not processing_window_data.missing_F_data, CALC_ERROR_CODE__MISSING_FDATA, "Missing F data"): 
                    break
                
                if not cycle_check( not processing_window_data.excluded_by_sample_number, CALC_ERROR_CODE__TOO_FEW_M_DATA, "Too few measurements"):
                    break

                                   
                assert hc_NED is not None
                
                # Check of balancing of the positions (to filter out partial measurements)
                NM_position_indeces, measurement_balanced =  du.determine_DIM_measurement_positions(dim_mag, hc_NED, vc_theo_rad, DIM_instr, processing_window_raw_data[dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION], processing_window_raw_data[dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION])
                NM_position_indeces_time_ordered = dict(sorted(NM_position_indeces.items(), key=lambda kv: min(kv[1], default=float('inf'))))
                
                processing_window_data.measurement_positions_order = list( NM_position_indeces_time_ordered.keys() ) 
                # NM_position_indeces = processing_window_data.NM_position_indeces
                # if not cycle_check( not processing_window_data.non_balalnce_measurement, CALC_ERROR_CODE__UNEQUAL_POSITIONS, "Non balanced measurement"):
                #    break
                
                if enable_cycle_plots:
    
                    # Creating a list of indeces, when DIM position changes (e.g. from Face-Left to Face-Right)
                    DIM_position_change_indeces = [] 
                    for pos_index in range(1,n):
                        if processing_window_raw_data[pos_index][dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION] != processing_window_raw_data[pos_index-1][dim_measurement_loader.COMPNAME__TELESCOPE_CARDINAL_DIRECTION] \
                        or processing_window_raw_data[pos_index][dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION] != processing_window_raw_data[pos_index-1][dim_measurement_loader.COMPNAME__TELESCOPE_FACE_POSITION] :
                            DIM_position_change_indeces += [ pos_index ]
                    
                    DIM_pos_index_separators = np.array(DIM_position_change_indeces) - 0.5
    
                    DIM_position_change_indeces_minus_1 = [x - 1 for x in DIM_position_change_indeces]                
                    DIM_pos_datetime_separators = (timelabel[DIM_position_change_indeces_minus_1] + timelabel[DIM_position_change_indeces] ) /2.0                        
                    
                    # mag_meas_point_datetimes = DIM_data_loader.M_readings[ mag_meas_file_records_indeces_all_calcvars ] [ dim_measurement_loader.COMPNAME__TIMESTAMP ]
                    mag_obs_dataset.plot_data_for_date_range( timelabel[0], timelabel[-1], mark_datetimes = timelabel, datetime_separators = DIM_pos_datetime_separators)
                    
    
                if  (calc_params.test_mode is not None) and (calc_params.test_mode == "measurements_after_30min"):
                    t_end = mu.exceltime_2_datetime(timelabel[-1])
                    if not cycle_check( t_end.minute > 31, CALC_ERROR_CODE__MISC_ERROR, "Not-after-30-min"): break
                
                
                if not cycle_check( measurement_balanced, CALC_ERROR_CODE__UNEQUAL_POSITIONS,  "Non balanced measurement" ): 
                    break
                    
                if cycle_skip_reason is not None: break
            
                add_result("duration",  (timelabel[-1] - timelabel[0])*1440) # Adding duration only for balanced measurements!

                assert isinstance(meas_method_info, list)
                meas_method_str = meas_method_info[1]

                # mag_sensor_index = calc_params.mag_sensor_index
                act_mag_sensor = DIM_instr.mag_sensors[calc_params.mag_sensor_index]

                starting_magfield   = copy.deepcopy(DIM_meas_cycle_magfield__common)
                starting_mag_sensor = copy.deepcopy(act_mag_sensor)

                solvers_variometer_I = solvers_variometer_D = None

                prelim_F,prelim_I,prelim_D = None,None,None
                solvers_H0 = solvers_D0 = solvers_Z0 = None
                # solvers_baseline_results_calc_approach  = None
                decl_incl_part_separate_results = {}

                variometer_mag_data, variomter_mag_data_with_window, variometer_mag_data_missing_record_indeces = DIM_data_loader.processing_cycle_parallel_mag_data_get(processing_window_data)
                   
                assert len(variometer_mag_data_missing_record_indeces) == 0, "Situation should be filtered out by missing_variometer_data"
                
                if not cycle_check( len(variometer_mag_data_missing_record_indeces) == 0, CALC_ERROR_CODE__MISSING_VARDATA, "No variometer data"): break

                assert not np.any( np.isnan(F_abs_pillar) ), "Situation should be filtered out by missing_F_data"

                if np.any(np.isnan(dim_mag)):
                    cu.warning("There is/are {} NAN values within {} DIM measurement!".format ( np.count_nonzero(np.isnan(dim_mag)), len(dim_mag)), True )
                    break
                

                # Results will be put in this table, and printed, using 'tabulate' library
                table = cu.tabulate_table_creator() # ["Param.name", "unit", "Value", "opt-min", "opt-max", "at-bound"]

                if (calc_params.varcorr_method!="none"):
                    # Will correct original ADS mag values with processing-window-time variation
                    # varcorr_mag_data = DIM_data_loader.processing_cycle_parallel_mag_data_get(processing_window_index, varcorr_data_parallel_mag_data_index, calc_params_index)
                    assert variometer_data_parallel_mag_data_index>=0

                    # if mu.structured_array_has_nans(varcorr_mag_data):
                    #    # It can happen, that we have an reference for only some of the data points.
                    #    # Can not correct one or more measuremen point due to missing correction data part
                    #    cu.warning("Variaton correction data is missing, skipping this cyce", True)
                    #    continue

                    if calc_params.varcorr_method != "LM/dyn":  # at 'LM/dyn', variation correction is done during LM cycles

                        dim_mag_varcorr_used = True
                        prev_calcvar = None if calc_params_index == 0 else calc_params_list[calc_params_index-1]

                        F_diff_check_limit = None
                        assert calc_params.pre_baseline_calculation_variometer_calibration_source in PRE_BASELINE_CALCULATION_VARIOMETER_CALIBRATION_SOURCES
                        assert calc_params.pre_baseline_calculation_dim_mag_sensor_parameters_source in PRE_BASELINE_CALCULATION_DIM_SENSOR_PARAM_SOURCES
                        # print("Using \"{}\" method for variation correction refinement...".format(calc_params.pre_baseline_calculation_variometer_calibration_source))

                        if calc_params.pre_baseline_calculation_dim_mag_sensor_parameters_source == "none":
                            varcorr_mag_sensor = copy.deepcopy(act_mag_sensor)
                        elif calc_params.pre_baseline_calculation_dim_mag_sensor_parameters_source == "calcvars-prev-results":
                            varcorr_mag_sensor = calc_params.last_results['mag_sensor']
                        elif calc_params.pre_baseline_calculation_dim_mag_sensor_parameters_source == "prev-calcvar-results":
                            assert prev_calcvar is not None
                            varcorr_mag_sensor = prev_calcvar.last_results['mag_sensor']
                        else:
                            assert False

                        # elif calc_params.pre_baseline_calculation_variometer_calibration_source == "GN":
                        #     assert False
                        #     varcorr_field = copy.deepcopy(DIM_meas_cycle_magfield__common)
                        #     mag_inverse_method_with_GN(dim_mag, hc_NED, vc_NED, F_abs_pillar, varcorr_mag_sensor, varcorr_field, b_plot_best_result_residuals = False )
                        #     del(varcorr_field)
                        # elif calc_params.pre_baseline_calculation_variometer_calibration_source in [ "LM-DI", "LM-D", "LM-|" ]:
                        #     assert False, "Not impldemented"
                        #     # mag_LM_inversions_for_DI_angles(dim_mag, hc_NED, vc_NED, F_abs_pillar, varcorr_mag_sensor, varcorr_field, b_plot_best_result_residuals = False )
                        if calc_params.pre_baseline_calculation_variometer_calibration_source == "user-defined":
                            assert mag_obs_dataset.convert_to_abs_data_def__preliminary is not None
                            prelim_F,prelim_I,prelim_D = gu.variometer_to_abs_data( calc_params.baseline_calc_approach, mag_obs_dataset, mag_obs_dataset.convert_to_abs_data_def__preliminary)
                            F_diff_check_limit = 100.0
                        elif calc_params.pre_baseline_calculation_variometer_calibration_source == "calcvars-prev-results":
                            # ToDo: in the first processing-window, baseline is not yet know, so the results below can be very false!
                            prelim_F,prelim_I,prelim_D = gu.variometer_to_abs_data( calc_params.baseline_calc_approach, mag_obs_dataset, calc_params_common, variometer_mag_data, calc_params.last_results["baseline"], allow_empty_baseline = True )
                            F_diff_check_limit = 60.0 # ToDo: levinni 1.5-re (vagy 6-result_array)
                        elif calc_params.pre_baseline_calculation_variometer_calibration_source == "prev-calcvars-results":
                            assert prev_calcvar is not None
                            prelim_F,prelim_I,prelim_D = gu.variometer_to_abs_data( prev_calcvar.baseline_calc_approach, mag_obs_dataset, variometer_mag_data, prev_calcvar.last_results["baseline"])
                            F_diff_check_limit = 1.2
                        else:
                            assert False

                        if prelim_F is not None:
                            F_diff_check = prelim_F - F_abs_pillar
                            max_F_diff_check = max(abs(F_diff_check))
                            if calc_params.do_self_checks:
                                assert abs(max_F_diff_check)<F_diff_check_limit, "F-diff: {}  limit:{} source:\"{}\"".format(max_F_diff_check, F_diff_check_limit, calc_params.pre_baseline_calculation_variometer_calibration_source)
                        if (prelim_I is None) or (prelim_D is None):
                            # At first processing-window, we accept not yet having full baseline system:
                            dim_mag_varcorr = np.zeros( len(hc_NED) )
                        else:
                            dim_mag_varcorr = du.calc_absmeas_variation_correction_values(calc_params.varcorr_method, hc_NED, vc_NED, varcorr_mag_sensor, F_abs_pillar, prelim_I, prelim_D )

                        if calc_params.varcorr_method != "none":
                            if (varcorr_prev_M_scale[calc_params_index] != None) and (abs(varcorr_mag_sensor.scale - varcorr_prev_M_scale[calc_params_index])  > 0.3):
                                cu.warning("Big change of M_scale between two cycles: {:.3f} => {:.3f}! The variation correction of the processing-window that do not (yet) have accurate M_scale may have error! (if M_scale sign is reversed, then correction is reversed also!)".format(varcorr_prev_M_scale[calc_params_index],act_mag_sensor.scale ), 30 )
                            varcorr_prev_M_scale[calc_params_index] = varcorr_mag_sensor.scale # This will be overwritten in GN processing-window
                            if np.any( np.isnan(dim_mag_varcorr)):
                                # Can not correct one or more measuremen point due to missing reference data
                                # Note: this situation probably handled by severalas lines above
                                assert False
                                cu.warning("Variaton correction data is not complete for this period, skipping this cyce")
                                continue

                            # ads_mag = ads_mag + dim_mag_varcorr

                        # Megj.: scale figyelembe van véve, ld. fentebb
                else:
                    dim_mag_varcorr = np.empty(n)
                    dim_mag_varcorr[:] = 0.0
                    dim_mag_varcorr_used = False

                assert not np.any(np.isnan(dim_mag))  # Checking 'dim_mag' after optional correction
                solver_params_at_the_bound_names = None # Some method returns the 'params_at_the_bound_names', that we also determine here later

                # By default 'ADS_cycle_magfield__D' and 'ADS_cycle_magfield__I' are pointers to the same, common 'DIM_meas_cycle_magfield__common'
                DIM_meas_cycle_magfield__common.indeces = list(range(n))
                ADS_cycle_magfield__Dpart = DIM_meas_cycle_magfield__common # Magfield for 'H' and 'D' modelling
                ADS_cycle_magfield__Ipart = DIM_meas_cycle_magfield__common # Magfield for 'Z' modelling
               

                if calc_params.method in [ "LM-DI", "LM-D", "LM-I", "LM-|" ]:
                    DIM_instr.ortho_error_arcsec = 0.0
                    assert calc_params.optimize_ortho_error == False
                    known_declination = not "D" in  calc_params.method[3:]
                    if known_declination:
                        DIM_meas_cycle_magfield__common.D = calc_params.known_D
                    known_inclination = not "I" in  calc_params.method[3:]
                    if known_inclination:
                        DIM_meas_cycle_magfield__common.I = calc_params.known_I
                    result_str, solver_params_at_the_bound_names, nof_outliers, resid_std, residuals = mag_LM_inversions_for_DI_angles(dim_mag, dim_mag_varcorr, hc_theo_rad, vc_theo_rad, F_abs_pillar, DIM_meas_cycle_magfield__common, calc_params, immediately_plot_1D_residual_for_each_cycle, hc_raw_tick, do_print_new_random_step_results, meas_method_str = meas_method_str, variometer_data = variometer_mag_data, known_inclination = known_inclination, known_declination = known_declination)
                    # Results are in 'DIM_meas_cycle_magfield__common' and in 'act_mag_sensor' (that is a link to  'calc_params.instr.mag_sensors[calc_params.mag_sensor_index]')
                    if calc_params.varcorr_method=="LM/dyn":
                        # calculating the final 'dim_mag_varcorr' now
                        # varcorr_mag_sensor = dim_instr.mag_sensor_model_class(M_offset = M_offset, M_scale = M_scale, epsilon = epsilon, delta = delta)
                        varcorr_mag_sensor = copy.deepcopy( DIM_instr.mag_sensors[calc_params.mag_sensor_index] )
                        dim_mag_varcorr = du.calc_absmeas_variation_correction_values(calc_params.varcorr_method, hc_NED, vc_NED, varcorr_mag_sensor, F_abs_pillar, variometer_mag_data[gu.COMPNAME_INCL_ANGLE], variometer_mag_data[gu.COMPNAME_DECL_ANGLE] )
                        dim_mag_varcorr_used = True

                    table.add('method', '', "{}", meas_method_str.replace("/","\n")  )

                elif calc_params.method[:6] == "LM-BL/":

                    # Levemberg-Marquardt method to get the baseline values directly from the optimization
                    # Note: it is already checked at the creation of the varcorr method instance
                    assert (calc_params.varcorr_method == "none"),"Variation correction should be not applied for this solver!"
                    
                    result_str, solver_params_at_the_bound_names, nof_outliers, resid_std, residuals, LM_results =  mag_LM_inversions_for_baseline(calc_params_common.baseline_calc_approach, has_refpos_data = False, sub_method = calc_params.method[6:], input_data = [ timelabel, dim_mag, hc_theo_rad, vc_theo_rad, hc_NED, vc_NED, F_abs_pillar], calc_var_index = calc_params_index, calc_var = calc_params, mag_obs_dataset = mag_obs_dataset, DIM_instr = DIM_instr,  b_plot_best_result_residuals = enable_cycle_plots, variometer_mag_data = variometer_mag_data,  b_plot__F_diff_vs_posID = False, show_results_table = enable_cycle_tables)

                    DIM_meas_cycle_magfield__common.I = np.mean(LM_results[0])
                    DIM_meas_cycle_magfield__common.D = np.mean(LM_results[1])
                    DIM_meas_cycle_magfield__common.F = np.mean(LM_results[2])
                    # LM_res_decbas = LM_results[3]

                    # This can be deleted later:
                    #if not calc_params.optimize_decbas:
                    #    if variometer_data_def.decbas is None:
                    #        assert LM_res_decbas == 0
                    #    else:
                    #        assert variometer_data_def.decbas == LM_res_decbas # 'ref_data_def.decbas' is used to recalculate H0/D0/Z0

                    solver_returned_A0 = LM_results[4]
                    solver_returned_B0 = LM_results[5]
                    solver_returned_Z0 = LM_results[6]

                    # table.add('method', '', "{}", meas_method_str.replace("/","\n")  )


                elif calc_params.method == "GN":
                    # Gauss-Newton method
                    residuals, varcorr_mag_sensor, dim_mag_varcorr  = mag_inverse_method_with_GN( dim_mag, dim_mag_varcorr, hc_NED, vc_NED, F_abs_pillar, DIM_instr, calc_params, DIM_meas_cycle_magfield__common, b_plot_best_result_residuals = True, variometer_mag_data = variometer_mag_data, mag_obs_dataset = mag_obs_dataset )
                    # "GN" has been executed more cycles, with gradual refinement of variation correction with the actual results.
                    varcorr_prev_M_scale[calc_params_index] = varcorr_mag_sensor.scale
                    resid_std = np.std(residuals)
                    table.add('method', '', "{}", meas_method_str.replace("/","\n"))
                    if immediately_plot_1D_residual_for_each_cycle:
                        plot_utils.plot_xy(range(0,len(residuals)), residuals, "means. point nr.", "resiudal (nT)", symbol='.', color='blue')

                elif calc_params.method in [ "RESID", "NM", "NM/old" ]:
                    # if not dim_measurement_loader.get_method_by_method_str(meas_method_str):
                    #     # Can not apply the analítical null method calculation for other types of measurements - skipping this cycle.
                    #     # 'ra_cycles[calc_cycle]' components should be filled by NANs at initialization. Just to be sure, we check it:
                    #     assert np.isnan( ra_cycles[calc_params_index].comp_array_get("cycle_nr")[processing_window_index] )
                    #     cu.warning('Skipping the null-method analytical calculation for cycle with method "{}"'.format(meas_method_str))
                    #     continue
                    assert (not calc_params.optimize_ortho_error) and (not calc_params.optimize_hc_eccentr_err) and (not calc_params.optimize_vc_eccentr_err)
                    assert calc_params.optimize_M_offset and calc_params.optimize_epsilon and calc_params.optimize_delta
                    assert not dim_mag_varcorr_used
                    errmsg, residual_method_results =  residual_method_anal_solver(calc_params_common.baseline_calc_approach, NM_position_indeces, dim_mag, hc_theo_rad, vc_theo_rad, hc_NED, vc_NED, F_abs_pillar, DIM_instr, act_mag_sensor, calc_params.test_mode, variometer_mag_data, mag_obs_dataset, calc_params_common, calc_params, meas_pos)
                    if errmsg is not None:
                        cycle_check(False, CALC_ERROR_CODE__SOLVERS_ERROR, "Residual-solver error: "+errmsg)
                        break
                    (solvers_variometer_I, solvers_variometer_D, I__I_part, I__D_part, D__I_part, D__D_part, Incl_series_indeces, Decl_series_indeces, F_Incl_avg, F_Decl_avg, resid_std, residuals, decl_incl_part_separate_results, solvers_H0,solvers_D0,solvers_Z0) = residual_method_results
                    
                    DIM_meas_cycle_magfield__common.I = np.mean(solvers_variometer_I)
                    DIM_meas_cycle_magfield__common.D = np.mean(solvers_variometer_D)

                    assert isinstance(DIM_meas_cycle_magfield__common.I, float)
                    if calc_params.method in [ "RESID", "NM"]:
                        # Separating D and I series modelled fields
                        # Commneted out in 2026.06.23, as it can result in problem if a cycle check fails, and process does not reach the recalculation line of thios F: 
                        # DIM_meas_cycle_magfield__common.F = None  # We clear this field just to ensure that we do not use 'F' of the common I/D magfield
                        DIM_meas_cycle_magfield__common.indeces = None   # To ensure that we do not use 'indeces' of the common I/D magfield
                        ADS_cycle_magfield__Dpart = copy.deepcopy(DIM_meas_cycle_magfield__common)
                        ADS_cycle_magfield__Ipart = copy.deepcopy(DIM_meas_cycle_magfield__common)
                        ADS_cycle_magfield__Dpart.F = F_Decl_avg
                        ADS_cycle_magfield__Dpart.D = D__D_part
                        ADS_cycle_magfield__Dpart.I = I__D_part  #  it was 'I__I_part' until 2025.12.15, and produced not accurate reults

                        ADS_cycle_magfield__Ipart.F = F_Incl_avg
                        ADS_cycle_magfield__Ipart.D = D__I_part #  it was 'D__D_part' until 2025.12.15, and produced not accurate reults
                        ADS_cycle_magfield__Ipart.I = I__I_part

                        ADS_cycle_magfield__Dpart.indeces= Decl_series_indeces
                        ADS_cycle_magfield__Ipart.indeces= Incl_series_indeces

                        datetime_Ipart = np.mean(timelabel[Incl_series_indeces])
                        datetime_Dpart = np.mean(timelabel[Decl_series_indeces])
                    else:
                        DIM_meas_cycle_magfield__common.F = np.mean(F_abs_pillar)
                        # 'ADS_cycle_magfield__D' and 'ADS_cycle_magfield__I' should be pointers to 'DIM_meas_cycle_magfield__common'
                        assert DIM_meas_cycle_magfield__common.I == ADS_cycle_magfield__Dpart.I
                        assert DIM_meas_cycle_magfield__common.D == ADS_cycle_magfield__Dpart.D

                    # epsilon_ID_diff_arcmin = angle_utils.rad_to_arcmin(epsilon_D-epsilon_I)

                    if immediately_plot_1D_residual_for_each_cycle:
                        plot_utils.plot_xy(range(0,len(residuals)), residuals, "means. point nr.", "resiudal (nT)", symbol='.', color='blue')

                else:
                    assert False, "Calculation method {} is not supported!".format(calc_params.method)


                # Ensuring that variation correction had been enough accurate:
                if calc_params.varcorr_method not in  [ "none", "LM/dyn" ]:
                    if abs(varcorr_mag_sensor.scale - act_mag_sensor.scale) > 0.05:
                        cu.warning('Variation correction time applied scale factor ({:.2f} nT/nT) and resulting scale factor ({:.2f} nT/nT) are quite different!'.format(varcorr_mag_sensor.scale, act_mag_sensor.scale), True)
                    if abs(varcorr_mag_sensor.epsilon_arcmin - act_mag_sensor.epsilon_arcmin) > 5:
                        cu.warning('Variation correction time applied epsilon  ({:.2f} arcmin) and resulting epsilon ({:.2f} arcmin) are quite different!'.format(varcorr_mag_sensor.epsilon_arcmin, act_mag_sensor.epsilon_arcmin), True)
                    if abs(varcorr_mag_sensor.delta_arcsec - act_mag_sensor.delta_arcsec) > 45:
                        cu.warning('Variation correction time applied delta ({:.2f} arcsec) and resulting delta ({:.2f} arcsec) are quite different!'.format(varcorr_mag_sensor.delta_arcsec, act_mag_sensor.delta_arcsec), True)

                calc_params.last_results['I'] = DIM_meas_cycle_magfield__common.I
                calc_params.last_results['D'] = DIM_meas_cycle_magfield__common.D
                calc_params.last_results['mag_sensor'] = copy.deepcopy(act_mag_sensor)

                calc_params.result_array[processing_window_index]["datetime-Incl"] = datetime_Ipart # Overwriting the common time
                calc_params.result_array[processing_window_index]["datetime-Decl"] = datetime_Dpart # Overwriting the common time
                
                # ------------------------
                # save_DI_opt_results(calc_params_index, meas_method_info[0], processing_window_index, cycle_name, cycle_file_id_info, timelabel )
                # def save_DI_opt_results():
                # ------------------------
                
                calc_params.result_array[processing_window_index]["error_code"] = CALC_ERROR_CODE__NONE  # Overwriting default error code
                

                add_result("azimark-hc_window-avg", angle_utils.user_angle(processing_window_data.azimark_readings_avg_angle ) )
                add_result("azimark-hc-change_mean-removed",     angle_utils.rad2arcsec(processing_window_data.azimark_readings_avg_angle) )

                # new_hc0_ned_SU = true_angle_calculator.h_theo_to_ned_from_new_SU(0.0)
                #new_hc0_ned_SD = true_angle_calculator.h_theo_to_ned_from_new_SD(0.0)
                # add_result_comp(True, True, "new-SU_hc(0)@NED", "deg",  math.degrees(new_hc0_ned_SU ), plot = ["new-SU_hc(0)@NED", "new-SD_hc(0)@NED", "new-SD_hc(0)@NED", "new-avg_hc(0)@NED"] )
                # add_result_comp(True, True, "new-SD_hc(0)@NED", "deg",  math.degrees(new_hc0_ned_SD ), plot = True )
                # add_result_comp(True, True, "new-avg_hc(0)@NED","deg",  math.degrees( (new_hc0_ned_SU+new_hc0_ned_SD)/2.0 ), plot = True )

                #if (true_angle_calculator.act_SU_hc_tick is None) or np.isnan(true_angle_calculator.act_SU_hc_tick):
                # add_result_comp(True, True, "new-SU",           "tick", true_angle_calculator.act_SU_hc_tick, plot = False)
                # add_result_comp(True, True, "new-SU",           "tick", 0.0 , plot = False) # true_angle_calculator.act_SU_hc_tick
                # add_result_comp(True, True, "new-SD",           "tick", 0.0, plot = False) # true_angle_calculator.act_SD_hc_tick

                # add_result_comp(True, True, "act-hc(0)@NED",    "deg",  math.degrees(true_angle_calculator.h_theo_to_ned(0.0) ), plot = True )

                
                add_result("azimark_WGS84",   math.degrees( processing_window_data.azimark_readings_true_calibrating_angle )  )
                # add_result("azimark_hc_theo", math.degrees( processing_window_data.azimark_readings_theo_angle )  )
                for hcvc in ["hc","vc"]:
                    add_result(f"azimark_{hcvc}-FL-FR-diff", angle_utils.rad2arcsec(processing_window_data.azimark_readings_avg_FL_FR_diff[hcvc]) )
                    
                add_result("vc_offset_from_azimark", math.degrees( processing_window_data.vc0_from_azimark_readings )   )

                # add_result_comp(True, True, "new-SU_vc", "deg",   0,0, plot = [ "new-SU_vc", "new-SD_vc-complement" ] ) # math.degrees(true_angle_calculator.act_SU_vc)
                # add_result_comp(True, True, "new-SD_vc-complement", "deg", 0.0, plot = True ) #   (360.0 - math.degrees(true_angle_calculator.act_SD_vc))

                # azi_SU_vc_tick = true_angle_calculator.act_SU_vc_tick if true_angle_calculator.act_SU_vc_tick is not None else np.nan
                # azi_SD_vc_tick_complement = DIM_instr.ENCODER_DIVISION - true_angle_calculator.act_SD_vc_tick + 2*dim_instr.rad_to_tick_raw(DIM_instr.vc0)  if true_angle_calculator.act_SD_vc_tick is not None else np.nan
                # add_result_comp(True, True, "new-SU_vc-tick", "tick",   azi_SU_vc_tick, plot = [ "new-SU_vc-tick", "new-SD_vc-complement-tick" ] )
                # add_result_comp(True, True, "new-SD_vc-complement-tick", "tick", azi_SD_vc_tick_complement  , plot = True )

                # new_SD_diff_tick = 0.0 # DIM_instr.vc_axis.rad_to_tick_raw  ( (math.tau - true_angle_calculator.act_SD_vc) - true_angle_calculator.act_SU_vc )

                # add_result_comp(True, True, "new-vc_diff", "tick",  new_SD_diff_tick, plot = True )
                
                if not cycle_check( resid_std < calc_params_common.resid_std_acceptance_limit, CALC_ERROR_CODE__TOO_BIG_RESID_STD, "Too big resdual error: {:.2f} nT".format(resid_std)):
                    break
                
                add_result("residual-stddev", resid_std)                  

                add_result("DIM-I",       math.degrees(DIM_meas_cycle_magfield__common.I), starting_val = math.degrees(starting_magfield.I), range_min = math.degrees(meas_pos.I_min) , range_max = math.degrees(meas_pos.I_max) )
                add_result("DIM-D",       math.degrees(DIM_meas_cycle_magfield__common.D), starting_val = math.degrees(starting_magfield.D), range_min = math.degrees(meas_pos.D_min),  range_max = math.degrees(meas_pos.D_max) )
                
                assert starting_magfield.F is not None # 2026.06.23: ADS_cycle_common.F can be deteled in a step above
                add_result("DIM-X",       ADS_cycle_magfield__Ipart.X(), starting_val = starting_magfield.X())
                add_result("DIM-Y",       ADS_cycle_magfield__Dpart.Y(), starting_val = starting_magfield.Y())
                add_result("DIM-Z",       ADS_cycle_magfield__Ipart.Z(), starting_val = starting_magfield.Z())

                # add_result("DIM-H",        ADS_cycle_magfield__Ipart.H())
                add_result("F",            (ADS_cycle_magfield__Dpart.F + ADS_cycle_magfield__Ipart.F)/2.0)
                add_result("M-offset",     act_mag_sensor.offset,         starting_mag_sensor.offset, act_mag_sensor.offset_min, act_mag_sensor.offset_max )
                add_result("M-scale",      act_mag_sensor.scale,          starting_mag_sensor.scale,  act_mag_sensor.scale_min, act_mag_sensor.scale_max )
                add_result("epsilon",      act_mag_sensor.epsilon_arcmin, starting_mag_sensor.epsilon_arcmin,   act_mag_sensor.epsilon_min_arcmin, act_mag_sensor.epsilon_max_arcmin )
                add_result("delta",        act_mag_sensor.delta_arcsec,   starting_mag_sensor.delta_min_arcsec, act_mag_sensor.delta_min_arcsec, act_mag_sensor.delta_max_arcsec )

                M_slope_from_variometer = np.empty( n )
                M_slope_from_variometer[:] = np.nan

                if variometer_mag_data is not None:
                    # Calculating residuals of the new measurement using the last baseline: it is used
                    # to track measurement errors
                    if calc_params.calc_residuals_from_last_result:
                        if calc_params.last_results["baseline"] is None:
                            residuals_from_last_result = np.empty( n )
                            residuals_from_last_result[:] = np.nan
                        else:
                            assert enable_cycle_plots, "This feature requires processing window plots to be enabled!"                            
                            prev_res_mag_sensor = calc_params.last_results["mag_sensor"]
                            abs_data_from_baseline_corrected_variometer_data = gu.variometer_to_abs_data( calc_params_common.baseline_calc_approach, mag_obs_dataset, calc_params_common, variometer_mag_data, calc_params.last_results["baseline"])
                            M_ideal_from_prev_resuls = amdp.dim_instrument_model(F_abs_pillar, vc_NED, hc_NED, abs_data_from_baseline_corrected_variometer_data.I, abs_data_from_baseline_corrected_variometer_data.D, prev_res_mag_sensor.epsilon(), prev_res_mag_sensor.delta(), prev_res_mag_sensor.offset, prev_res_mag_sensor.scale)
                            residuals_from_last_result = dim_mag - M_ideal_from_prev_resuls
                    

                    calc_params.last_results["baseline"] = calc_baseline_values(calc_params_common.baseline_calc_approach,  mag_obs_dataset, variometer_mag_data, F_abs_pillar, ADS_cycle_magfield__Dpart, ADS_cycle_magfield__Ipart, calc_params_common)
                    baseline_values = calc_params.last_results["baseline"]
                    
                    def post_check_error(err_code, err_msg):
                        cu.warning(err_msg)
                        if calc_params.result_array[processing_window_index]["error_code"] == CALC_ERROR_CODE__NONE:
                            calc_params.result_array[processing_window_index]["error_code"] = err_code
                        # ToDo: consider invalidating I/D result at the end 
                    
                    # if (solvers_baseline_results_calc_approach is not None) and (solvers_baseline_results_calc_approach == calc_params.baseline_calc_approach):
                    if True:
                        if solvers_H0 is not None:
                            if not cycle_check(abs(solvers_H0-baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR])<0.1 * 100, CALC_ERROR_CODE__POST_CHECK_H0, f"H0 check error: {solvers_H0} {baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR]}"): break
                        if (calc_params_common.baseline_calc_approach == "D0-angle" ):
                            diff_between_D0_arcsec = abs( math.degrees(solvers_D0-baseline_values[gu.COMPNAME__BASELINE__D_ANGLE])*3600)
                            if calc_params.do_self_checks:
                                assert diff_between_D0_arcsec < 20, diff_between_D0_arcsec # ToDo: set back to 5
                                # ToDo: call 'post_check_error()' instead
                        if solvers_Z0 is not None:
                            if calc_params.do_self_checks:
                                assert abs(baseline_values[gu.COMPNAME__BASELINE__Z_SENSOR]-solvers_Z0) < 0.01
                                # if abs (Baseline_Z - solver_returned_Z0) > 0.1: cu.warning("Difference in different calculated Z baselines: solver:{:.1f} / post-calc:{:.1f}".format(solver_returned_Z0, Baseline_Z), True)
                                # ToDo: call 'post_check_error()' instead
                                
                    if calc_params.baseline_tempco_corr is not None:

                        correcting_temp = aux_data[aux_data_idx_for_this_cycle][calc_params.baseline_tempco_corr_source]

                        assert len(calc_params.baseline_tempco_corr['H0']) == 2
                        M = calc_params.baseline_tempco_corr['H0'][0]
                        A = calc_params.baseline_tempco_corr['H0'][1]
                        assert calc_params.baseline_tempco_corr_source in aux_data.dtype.names

                        baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR] -= correcting_temp*M - A

                        M = calc_params.baseline_tempco_corr['D0_min'][0]
                        A = calc_params.baseline_tempco_corr['D0_min'][1]

                        baseline_values[gu.COMPNAME__BASELINE__D_ANGLE] -= correcting_temp*M - A

                        M = calc_params.baseline_tempco_corr['Z0'][0]
                        A = calc_params.baseline_tempco_corr['Z0'][1]

                        baseline_values[gu.COMPNAME__BASELINE__Z_SENSOR] -= correcting_temp*M - A
                        
                    # --------------------------------------------------------------------
                    # Calculating full field from recently calibrated variometer
                    abs_data_from_baseline_corrected_variometer_data = gu.variometer_to_abs_data( calc_params_common.baseline_calc_approach, mag_obs_dataset, calc_params_common, variometer_mag_data, baseline_values)


                    if calc_params_common.baseline_calc_approach == "D0-angle":
                        assert XH_compname in ['X','H']
                        assert YD_compname in ['D','Y']
                        add_result("H0", baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR] )
                        
                        D0_baseline_arcmin = np.degrees(baseline_values[gu.COMPNAME__BASELINE__D_ANGLE])*60.0
                        add_result("D0-angle", D0_baseline_arcmin) 
                        
                        # H_baseline_corrected = abs_data_from_baseline_corrected_variometer_data.H()
                        # H_baseline_corrected = 13200 # qqq
                        
                        # # D0_nT_measured           =  H_baseline_corrected * np.tan(  abs_data_from_baseline_corrected_variometer_data.D - baseline_values[gu.COMPNAME__BASELINE__D_ANGLE]  ) 
                        # # D0_nT_baseline_corrected =  H_baseline_corrected * np.tan(  abs_data_from_baseline_corrected_variometer_data.D )
                        
                        # D0_nT = np.mean(H_baseline_corrected) * D0_baseline_arcmin / 3437.75
                        # # 3437.75 being the number of arcminutes per radian
                        # D0_nT = np.mean(D0_nT)
                        # add_result("D0_nT_from_angle",  np)

                    elif calc_params_common.baseline_calc_approach == "Y0-nT":
                        assert XH_compname in ['X','H']
                        assert YD_compname in ['D','Y']
                        add_result(XH_compname + "0",  baseline_values[gu.COMPNAME__BASELINE__HX_SENSOR] )
                        add_result("D0-angle",         np.nan )
                        add_result(YD_compname + "0",  baseline_values[gu.COMPNAME__BASELINE__DY_SENSOR] )
                    
                    elif calc_params.baseline_calc_approach == "dIdD":  
                        assert mag_obs_dataset.instr_frame() == "dIdD"
                        
                        
                    else:
                        assert False, calc_params.baseline_calc_approach

                    add_result("Z0",  baseline_values[gu.COMPNAME__BASELINE__Z_SENSOR] )
                    
                    # Removed: 2026.01.27
                    # if not calcvar_has_valid_cycle_result[calc_params_index]:
                    #     mag_comps = [ mag_obs_dataset.instr_frame()[0]+'0'  ]  # H0 or X0
                    #     concatenated_plot_y_components += [  ]  # ???
                    #     if calc_params.baseline_calc_approach == "D0-angle":
                    #         mag_comps +=  [   "D0_nT_from_angle", "D0-angle"   ]
                    #     else:
                    #         mag_comps += [ YD_compname+'0' ]
                    #     mag_comps += [ 'Z0' ]
                    #     concatenated_plot_y_components += [ mag_comps ]



                    # if variometer_data_def.frame == 'XYZ':
                    #     add_result_comp(True, True, "X0",  "nT",  Baseline_A + (0  if (variometer_data_def.X0_applied_at_file_reading is None) or (not variometer_data_def.readd_baselines_applied_at_file_reading)  else  variometer_data_def.X0_applied_at_file_reading)  )
                    #     add_result_comp(True, True, "Y0",  "nT",  Baseline_B + (0  if (variometer_data_def.Y0_applied_at_file_reading is None) or (not variometer_data_def.readd_baselines_applied_at_file_reading)  else  variometer_data_def.Y0_applied_at_file_reading)  )
                    # elif variometer_data_def.frame == 'HDZ':  # HDZ oriented variometer, with known 'decbas'
                    #     add_result_comp(True, True, "H0",  "nT",  Baseline_A + (0 if (variometer_data_def.H0_nT_applied_at_file_reading is None) or (not variometer_data_def.readd_baselines_applied_at_file_reading)  else variometer_data_def.H0_nT_applied_at_file_reading) )
                    #     add_result_comp(True, True, "D0",  "nT",  Baseline_B + (0 if (variometer_data_def.D0_nT_applied_at_file_reading is None) or (not variometer_data_def.readd_baselines_applied_at_file_reading)  else variometer_data_def.D0_nT_applied_at_file_reading) )
                    #     if calc_params.method[:6] == "LM-BL/":
                    #         # Checking if directly cacluatled baselines are fit to baselines determined by the indierct way:
                    #             if abs (Baseline_A - solver_returned_A0) > 0.1: cu.warning("Difference in different calculated X/H baselines: solver:{:.1f} / post-calc:{:.1f}".format(solver_returned_A0, Baseline_A), True)
                    #             if abs (Baseline_B - solver_returned_B0) > 0.1: cu.warning("Difference in different calculated Y/D baselines: solver:{:.1f} / post-calc:{:.1f}".format(solver_returned_B0, Baseline_B), True)
                    # else:
                    #     assert False

                        
                    # add_result_comp(True, True, "Fdiff-avg","nT",     DIM_meas_cycle_magfield__common.F  -  np.mean( variometer_data_def.calc_F(variometer_mag_data, [Baseline_A, Baseline_B, Baseline_Z] ) ), plot = True  )
                    # add_result_comp(True, True, "Fdiff-stddev","nT",  np.std( DIM_meas_cycle_magfield__common.F   -  variometer_data_def.calc_F(variometer_mag_data, [Baseline_A, Baseline_B, Baseline_Z] ) ), plot = True  )

                    ref_mag_data_avg_TotalH_at_Incl_series = np.mean( gu.XY_2_H(variometer_mag_data[ YD_compname ][ADS_cycle_magfield__Ipart.indeces], variometer_mag_data[ YD_compname ][ADS_cycle_magfield__Ipart.indeces]))
                    #add_result_comp(True, True, "Horiz-diff",  "nT", ADS_cycle_magfield__Ipart.H() - ref_mag_data_avg_TotalH_at_Incl_series, plot = True )

                    #if gu.COMPNAME_INCL_ANGLE in variometer_mag_data.dtype.names:
                    #    add_result_comp(True, True, "Incl-diff", '"',  math.degrees(ADS_cycle_magfield__Ipart.I - np.mean(variometer_mag_data[gu.COMPNAME_INCL_ANGLE][ADS_cycle_magfield__Ipart.indeces]))*3600, plot = True )  # result_tabulate_data += [ [ "I-diff", "arcsec", "{:.4f}".format(math.degrees(DIM_meas_cycle_magfield__common.I - ref_I_cycle_avg)*3600)  ] ]
                    #if gu.COMPNAME_DECL_ANGLE in variometer_mag_data.dtype.names:
                    #    add_result_comp(True, True, "Decl-diff", '"',  math.degrees(ADS_cycle_magfield__Dpart.D - np.mean(variometer_mag_data[gu.COMPNAME_DECL_ANGLE][ADS_cycle_magfield__Ipart.indeces]))*3600, plot = True )  #


                    #if 'T1' in variometer_mag_data.dtype.names:
                    #    add_result_comp(True, True, "Variometer-T1", 'K',  np.mean(variometer_mag_data['T1']), plot = True )
                    #    add_result_comp(True, True, "Variometer-T2", 'K',  np.mean(variometer_mag_data['T2']), plot = True )
                    
                    # If we know discrete I/D values for meas. point, then we use that values (solvers_variometer_I / solvers_variometer_D)
                   
                    if (solvers_variometer_I is not None) and calc_params.do_self_checks:
                        # 'solvers_variometer_I' and  'solvers_variometer_D' is returned by the solver -> we should get the same here
                        max_I_diff_arcsec = 3600. *  math.degrees(max(abs(abs_data_from_baseline_corrected_variometer_data.I - solvers_variometer_I))) 
                        if  (max_I_diff_arcsec > 13.0):  
                            post_check_error(CALC_ERROR_CODE__POST_CHECK_I_1, "Big I difference: {:.2f} arcsec".format(max_I_diff_arcsec) )  
                    if (solvers_variometer_D is not None) and calc_params.do_self_checks:
                        max_D_diff_arcsec = 3600. * math.degrees( max(abs(abs_data_from_baseline_corrected_variometer_data.D - solvers_variometer_D)) )
                        if (max_D_diff_arcsec > 15):  
                            post_check_error(CALC_ERROR_CODE__POST_CHECK_I_2, "I check error 1: {:.2f} arcsec".format(max_D_diff_arcsec) )  
                        
                        
                    if calc_params_common.variometer_slopecalc_halfwindow_sec  != 0:
                        # We have variometer recording around each DIM measurement point's time
                        # We convert these variometer data to the DIM measurement positions
                        
                        N_window = 2 * calc_params_common.variometer_slopecalc_halfwindow_sec + 1
                        M_variometer_at_dim = np.empty( (n,N_window) ) # QQQ: lehetne egy dimenziósítani!
                        variometer_in_window_F,variometer_in_window_I,variometer_in_window_D = gu.variometer_to_abs_data( calc_params.baseline_calc_approach, mag_obs_dataset, variomter_mag_data_with_window, calc_params.last_results["baseline"])
                        LSQ_X__seconds = np.linspace(-dim_measurement_loader.qqq_variometer_data_halfwindows_sec, +dim_measurement_loader.qqq_variometer_data_halfwindows_sec, N_window)
                        for iwindow_sample in range(N_window):
                           M_variometer_at_dim[:,iwindow_sample] = amdp.dim_instrument_model(F_abs_pillar, vc_NED, hc_NED, variometer_in_window_I[:,iwindow_sample], variometer_in_window_D[:,iwindow_sample], act_mag_sensor.epsilon(), act_mag_sensor.delta(), act_mag_sensor.offset, act_mag_sensor.scale)
                        M_diff_check = M_variometer_at_dim[:,calc_params_common.variometer_slopecalc_halfwindow_sec]-dim_mag  # Difference between the center sample from the window (sample taken at DIM meas time) and DIM mag
                        if np.mean(M_diff_check) > 20.0:
                            post_check_error(CALC_ERROR_CODE__POST_CHECK_1,  "M_diff_mean:{}  M-diff:{}".format(np.mean(M_diff_check), M_diff_check) )

                    
                        # Calculating field change (slope) at each measurement point, using short term variation at measurement times
                        # In addition, checking result by comparing variometer average value (offset) and DIM value 
                        # M_variometer_at_dim/variometer_in_window_I/variometer_in_window_D 1st dimension: cycle point index
                        #                                                                   2nd dimension: index within window
                        
                        error_reported = False
                        for i in range(n):
                            if np.count_nonzero(~np.isnan(M_variometer_at_dim[i])) > 2:  # Need at least 3 variometer data
                                offset, M_slope_from_variometer[i], resid_info = mu.LSQ( LSQ_X__seconds,  M_variometer_at_dim[i], "o" )
                                # plt.plot(LSQ_X__seconds,  M_variometer_at_dim[i])
                                # plt.plot(LSQ_X__seconds,  offset + M_slope_from_variometer[i] * LSQ_X__seconds, "-", color = "red" )
                                # plt.show()
                                # cu.wait_for_enter()
                                if calc_params.do_self_checks: 
                                    if abs(offset - dim_mag[i]) > 150.0:
                                        if not error_reported: 
                                            # Probably there is some disturbance in DIM or in variometer data
                                            post_check_error(CALC_ERROR_CODE__POST_CHECK_2, "Big difference between variometer ({:.1f} nT) data and DIM ({:.1f} nT) result".format(offset,dim_mag[i]) )  
                                            error_reported = True
                            else:
                                # We can be here when having non continuous variometer and/or scalar data (e.g. F sampled at 5 Hz, and not interpolated to 1 Hz)
                                assert np.isnan(M_slope_from_variometer[i])
                                
                    if calc_params.calculated_null_positions_file is not None:
                        data_for_time_reduction = { "hc_NED":hc_NED, "vc_NED":vc_NED, "mag-sensor": act_mag_sensor, "magfield": abs_data_from_baseline_corrected_variometer_data }
                        # Comparing of results for theo and NED coordinates can be interesting
                        rotated_angles_theo, fixed_angles_theo, scale_check_measurement_theo,   zeropos_LSQ_resid_std_theo   = calc_null_method_positions(NM_position_indeces_time_ordered, dim_mag, hc_theo_deg, vc_theo_deg, data_for_time_reduction)
                        rotated_angles_NED,  fixed_angles_NED,  scale_check_measurement_NED,    zeropos_LSQ_resid_std_NED    = calc_null_method_positions(NM_position_indeces_time_ordered, dim_mag, np.degrees(hc_NED), np.degrees(vc_NED), data_for_time_reduction)
                        calc_params.calculated_null_positions[processing_window_index] = [ rotated_angles_theo, fixed_angles_theo, scale_check_measurement_theo, zeropos_LSQ_resid_std_theo, rotated_angles_NED ]


                for decl_incl_part_separate_result in decl_incl_part_separate_results:
                     # Extra results in case on NM solver: M_offset_D / M_offset_I / epsilon_D / epsilon_I
                #     add_result(True, True, decl_incl_part_separate_result["name"]+'-D', decl_incl_part_separate_result["unit"], decl_incl_part_separate_result["value-D"], optimized_param = True, starting_val = None, range_min = decl_incl_part_separate_result["range_min"], range_max = decl_incl_part_separate_result["range_max"], plot = True)
                #     add_result(True, True, decl_incl_part_separate_result["name"]+'-I', decl_incl_part_separate_result["unit"], decl_incl_part_separate_result["value-I"], optimized_param = True, starting_val = None, range_min = decl_incl_part_separate_result["range_min"], range_max = decl_incl_part_separate_result["range_max"], plot = True)
                      add_result(decl_incl_part_separate_result["name"]+'_D/I-diff', decl_incl_part_separate_result["value-D"] - decl_incl_part_separate_result["value-I"])
                #     if not calcvar_has_valid_cycle_result[calc_params_index]:
                #         concatenated_plot_y_components += [ [decl_incl_part_separate_result["name"]+'-D',decl_incl_part_separate_result["name"]+'-I' ]  ]

                # if calc_params.optimize_hc_eccentr_err:
                #     add_result_comp(True,True,  "hc-ecc_phase",   "°",  math.degrees(DIM_instr.hc.eccentr_error_phase), True, plot = True )
                #     add_result_comp(True,True,  "hc-ecc_ampl",   '"',  angle_utils.rad_to_arcsec(DIM_instr.hc.eccentr_error_ampl), True, range_min=0.0, range_max = angle_utils.rad_to_arcsec(DIM_instr.hc.max_eccentr_error_ampl), plot = True )
                # if calc_params.optimize_vc_eccentr_err:
                #     add_result_comp(True,True,  "vc-ecc_phase",   "°",  math.degrees(DIM_instr.vc.eccentr_error_phase), True, plot = True )
                #     add_result_comp(True,True,  "vc-ecc_ampl",   '"',  angle_utils.rad_to_arcsec(DIM_instr.vc.eccentr_error_ampl), True, range_min=0.0, range_max = angle_utils.rad_to_arcsec(DIM_instr.vc.max_eccentr_error_ampl), plot = True )


                # if aux_data is not None:
                #     for aux_comp_name in aux_data.dtype.names[1:]:
                #         assert False, "Use average value for the entire cycle"
                #         aux_value = aux_data[aux_data_idx_for_this_cycle][aux_comp_name]
                #         add_result_comp(False,True,  aux_comp_name,   '',  aux_value , False, plot = True, comp_prec = 2 )

                # add_result_comp(True, True, "refmark-WGS84",     "gon",  angle_utils.rad2gon(true_angle_calculator.last_measured_refpos_hc_WGS84[ DIM_instr.hc_angle_reference_sensor_pos()])  )
                
                # fourhour_index= int(hour_of_day) // 4
                # add_result_comp(False, True, "fourhours-of-day",    "_",  fourhour_index, comp_prec = 0)    
                
                
                # Storing (and saving) point-wise data
                # TdDo: reise this part
                NMposition_mean_removed_residuals = du.NMposition_mean_removed_residuals_calc(residuals) # For 'M_mean_removed_resid_nT' component
                
                ra_measpoint_indeces = slice(calc_params.nr_of_ra_measpoints_records, calc_params.nr_of_ra_measpoints_records+n)
                calc_params.ra_measpoints["datetime"][ra_measpoint_indeces] = timelabel
                calc_params.ra_measpoints["win_nr"][ra_measpoint_indeces]  = processing_window_index
                calc_params.ra_measpoints["pos_ID"][ra_measpoint_indeces]  = processing_window_raw_data[dim_measurement_loader.COMPNAME__POS_ID]
                calc_params.ra_measpoints["win_sample_nr"][ra_measpoint_indeces]  = np.arange(n)
                calc_params.ra_measpoints["hc_THEO"][ra_measpoint_indeces] = np.degrees(hc_theo_rad)
                calc_params.ra_measpoints["vc_THEO"][ra_measpoint_indeces] = np.degrees(vc_theo_rad)
                calc_params.ra_measpoints["hc_NED"][ra_measpoint_indeces]  = np.degrees(hc_NED)
                calc_params.ra_measpoints["vc_NED"][ra_measpoint_indeces]  = np.degrees(vc_NED)
                calc_params.ra_measpoints["M"][ra_measpoint_indeces]       = dim_mag
                calc_params.ra_measpoints["M_raw_resid_nT"][ra_measpoint_indeces]  =  residuals
                calc_params.ra_measpoints["M_mean_removed_resid_nT"][ra_measpoint_indeces]  = NMposition_mean_removed_residuals
                
                if dim_measurement_loader.COMPNAME__MAG_NOISE in  processing_window_raw_data.dtype.names:
                    calc_params.ra_measpoints["M_meas_noise"][ra_measpoint_indeces] = processing_window_raw_data[dim_measurement_loader.COMPNAME__MAG_NOISE]                    
                else: 
                    calc_params.ra_measpoints["M_meas_noise"][ra_measpoint_indeces] = np.nan                    

                if dim_measurement_loader.COMPNAME__MAG_SLOPE in processing_window_raw_data.dtype.names:
                    calc_params.ra_measpoints["M_meas_slope"][ra_measpoint_indeces] = processing_window_raw_data[dim_measurement_loader.COMPNAME__MAG_SLOPE]
                else:
                    calc_params.ra_measpoints["M_meas_slope"][ra_measpoint_indeces] = np.nan
                if  do_tick_mod_4_analysis: 
                    calc_params.ra_measpoints["hc_tick_mod_4"][ra_measpoint_indeces] = hc_raw_tick % 4
                    calc_params.ra_measpoints["hc_tick_mod_4"][ra_measpoint_indeces] = vc_raw_tick % 4
                    
                    
                # calc_params.ra_measpoints["M_slope_from_variometer"][ra_measpoint_indeces]  = M_slope_from_variometer
                # calc_params.ra_measpoints["F"][ra_measpoint_indeces]  = residuals
                #for i in range(n):
                    #if solvers_variometer_I is not None:
                    #    measpoint_I, measpoint_D =  solvers_variometer_I[i], solvers_variometer_D[i]
                    #else:
                    #    measpoint_I = DIM_meas_cycle_magfield__common.I
                    #    measpoint_D = DIM_meas_cycle_magfield__common.D
                    # epsilon = act_mag_sensor.epsilon()
                    # delta = act_mag_sensor.delta()
                    # offset = act_mag_sensor.offset, 
                    # scale = act_mag_sensor.scale, 
                    # dim_mag_varcorr[i], 
                    # variometer_mag_data[XH_compname][i], variometer_mag_data[YD_compname][i], variometer_mag_data["Z"][i]  )
                    # F_abs_pillar

                calc_params.nr_of_ra_measpoints_records  +=  n

                # avg-sum of point wise residuals (for each pos-id):
                assert n == len(processing_window_raw_data)
                for i in range(n):
                    pos_id = processing_window_raw_data[dim_measurement_loader.COMPNAME__POS_ID][i]
                    if not pos_id in calc_params.resid_per_pos_id_avgsum:
                        calc_params.resid_per_pos_id_avgsum[pos_id] = residuals[i]
                        calc_params.resid_per_pos_id_count[pos_id] = 1
                    else:
                        calc_params.resid_per_pos_id_avgsum[pos_id] += residuals[i]
                        calc_params.resid_per_pos_id_count[pos_id]  += 1

                # Printing out results in a tabulated form:
                if enable_cycle_tables:
                    print(tabulate.tabulate(table.result_tabulate_data, headers=table.headers, tablefmt="pretty"))  # tablefmt="simple_outline"

                params_at_the_bound_names = table.params_at_the_bound_names
                if (solver_params_at_the_bound_names is not None) and (len(params_at_the_bound_names) != len(solver_params_at_the_bound_names)):
                    # 'LM-BASELINE' can have params at bound for H0/D0/Z0, that is missing in table.params_at_the_bound_names
                    for param_at_bound_name in solver_params_at_the_bound_names:
                        if not param_at_bound_name in params_at_the_bound_names:
                            params_at_the_bound_names += [ "SOLVER[{}]".format(param_at_bound_name) ]

                    # assert False, "{} <> {}".format(params_at_the_bound_names,check_params_at_the_bound_names)
                if len(params_at_the_bound_names) > 0:
                    cu.warning('Some optimized parameters are at the border of optimization range: "{}"'.format(params_at_the_bound_names), 30)

                #if ( len(hc_raw_tick) > len(residuals) ):
                #    cu.warning("One or more outliers have been removed from calcualation - not using this cycle for mod 4 plot!")

                if enable_cycle_plots:
                    
                    
                    fig = plt.figure(figsize=(preferences.FIGURE_SIZE[0] * 0.7, preferences.FIGURE_SIZE[1] * 2 * 0.7) , layout = ['constrained', 'compressed', 'tight', 'none', None][2])
                    fig.suptitle('Residuals')
                    plt.subplot(2,1,1)
                    if dim_mag_varcorr_used:
                        plot_utils.plot_xy( range(n),  [ -dim_mag_varcorr,residuals]  , "sample count", [ "Calculated correction values (nT)", "Resid (nT)"], ['x', 'o'], color=['blue', 'red'], close_plot = False )
                        plt.subplot(2,1,2)
                        plot_utils.plot_datetime_y( timelabel,  [ -dim_mag_varcorr, residuals]  ,"nT" , label = [ "Calculated correction values (nT)", "Resid (nT)"], color=['blue', 'red' ], fmt=['x', 'o'] )
                    else:
                        Y = [residuals]
                        colors=[ 'red']
                        symbols = ['o']
                        series_names = None
                        if calc_params.calc_residuals_from_last_result:
                            Y += [residuals_from_last_result]
                            colors += ['blue']
                            symbols += ['x']
                            series_names = [ "Resid. from new res.", "Resid. from prev. res" ]
                        
                        plot_utils.plot_xy( range(n),  Y , "sample count", "Residuals (nT)", symbols, colors, close_plot = False, vertical_separators = DIM_pos_index_separators, label = series_names )
                        plt.subplot(2,1,2)
                        
                        datetimeformatter = "%H:%M:%S"
                        plot_utils.plot_datetime_y( timelabel,  [ residuals ]  ,"nT" , label = ["Resid (nT)"], color=[ 'red' ], fmt=['o'], datetimeformatter = datetimeformatter, vertical_separators = DIM_pos_datetime_separators )

                    plt.show()

                # 2-D plot
                if do_2D_residual_plot_for_each_cycle :
                    plot_utils.contour_plot(hc_theo_deg, vc_theo_deg, residuals)
                    # hf = plt.figure()
                    # ha = hf.add_subplot(111, projection='3d')
                    # ha.plot_trisurf(hc_deg, vc_deg, residuals)
                    # plt.xlabel('hc (°)')
                    # plt.ylabel('vc (°)')
                    # plt.show()

                DIM_meas_cycle_magfield__common.F = (ADS_cycle_magfield__Dpart.F + ADS_cycle_magfield__Ipart.F) / 2.0

            # End of dummy single pass cycle - we get here in case of errors (handled by 'break'), and in case of results
            if cycle_skip_reason is not None:
                # E.g. solver error or no reference mark
                print("  Not calculated due to problem: "+cycle_skip_reason)
                # cu.warning(cycle_skip_reason)                
           
            if (enable_cycle_tables or enable_cycle_plots) and (continuous_mode=="waitforenter"):
                resp = input("[Enter]=continue, b=break, a=abort, c=run-continuously: ") # ToDo: attenni du.user_input-result_array
                if resp == 'b':
                    continuous_mode="break"
                    # warning("A program további része nincs felkeszulve arra, hogy ido előtt megálljon a számítás, ugyanis az összefűzött ciklusok összes eredményére le lett foglalva a hely")
                if resp == 'a': sys.exit(1)
                if resp == 'c': continuous_mode="cont"

            if continuous_mode=="break":
                break

        # End of calc. variation loop
        if continuous_mode=="break":
            break

    # End of cycle loops
    # ----------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------

    if enable_post_plots or (len(calc_params_common.post_processing_steps)>0):
        cu.print_header("End of processing windows analysis, saving/printing/plotting results together")
    else:
        print("--- done ---")

    
    nr_of_all_caclvars_valid_results = 0
    min_result_time = 99999999
    max_result_time = 0
    
    

    for icalcvar, calc_params in enumerate(calc_params_list):
        
        calc_params.nr_of_valid_results =  np.count_nonzero( calc_params.result_array["error_code"] == CALC_ERROR_CODE__NONE  ) 
        if calc_params.nr_of_valid_results == 0:
            print("Calc params \"{}\" has any result!".format( calc_params.name))
            # print("There aren't any results in {:} calculation parameter set", icalcvar + 1)
        else:
            nr_of_all_caclvars_valid_results += calc_params.nr_of_valid_results
            if calc_params_common.baseline_calc_approach == "D0-angle":
                H = np.sqrt( np.square(calc_params.result_array["DIM-X"] ) + np.square(calc_params.result_array["DIM-Y"] ) )
                H_avg = np.nanmean( H )
                # Simple conversion: calc_params.result_array["D0_nT_from_angle"] = H_avg * calc_params.result_array["D0-angle"] / 3437.75 # 3437.75 being the number of arcminutes per radian
                calc_params.result_array["D0_nT_from_angle"] = H_avg * np.tan( angle_utils.arcmin_to_rad(calc_params.result_array["D0-angle"]) ) 
                
       
        # Removing means
        for comp in remove_mean_comps:
            if comp in calc_params.result_array.dtype.names:
                if not np.all( np.isnan(calc_params.result_array[comp])) :
                    mean = np.nanmean(calc_params.result_array[comp])
                    calc_params.result_array[comp][:] = calc_params.result_array[comp][:] - mean
            else:
                if calc_params_list[icalcvar].nr_of_valid_results>0:
                    assert False, f'Component "{comp}" is marked to be mean removed, but it is not within results'
                    

        # Saving results:
            
        result_array_comps_params["win_nr"].format_mask = "%{}u".format( len(str(DIM_data_loader.number_of_processing_windows))  ) 
        result_array_comps_params["input_file_win_idx"].format_mask = "%{}u".format( len(str( np.max(calc_params.result_array["input_file_win_idx"]) ))  ) 
        
        highest_stddev = np.nanmax( calc_params.result_array["residual-stddev"]) 
        highest_stddev_int_part_len = 3 if np.isnan(highest_stddev) else len( str(  int(highest_stddev)   )    )
        stddev_decimal_digits = 2
        result_array_comps_params["residual-stddev"].format_mask = "%{}.{}f".format(highest_stddev_int_part_len + 1 + stddev_decimal_digits, stddev_decimal_digits)
        
        min_result_time = min( np.min(calc_params.result_array["datetime-start"]), min_result_time)
        max_result_time = max( np.max(calc_params.result_array["datetime-end"]),   max_result_time)
        mu.savetxt(calc_params.result_file_name, calc_params.result_array, result_array_comps_params)
        saved_files += [calc_params.result_file_name]
        

        # Removing last unused records from 'ra_measpoints'
        calc_params.ra_measpoints = calc_params.ra_measpoints[:calc_params.nr_of_ra_measpoints_records]
                
        if calc_params.measpoints_result_file_name is not None:
            
            # For nice header/field format, see: https://gist.github.com/Jhsmit/3b066e1cecd27f769ee679094f7e14cd
            np.savetxt( calc_params.measpoints_result_file_name,calc_params.ra_measpoints, header = ra_measpoints_header, comments = "", fmt=ra_measpoints_fmt  )
            saved_files += [calc_params.measpoints_result_file_name]
            

        # Saving null method calculated positions:            
        if calc_params.calculated_null_positions_file is not None:
           with open(calc_params.calculated_null_positions_file, "w") as f: 
               
               # vc_180deg__differences = { "FL": [], "FR": [] }
               # for window_index in calc_params.calculated_null_positions.keys():  # processing window indeces are in the keys
               #     window = calc_params.processing_window_data[window_index]
               #     f.write("-- Null method positions of calculation window {} - date-range: {}...{} --\n".format( window.name, window.first_M_timelabel, window.last_M_timelabel ) )
               #     rotated_angles, fixed_angles, scale_check_data, zeropos_LSQ_residstd = calc_params.calculated_null_positions[window_index]
               #     for pos, rotated_angle in rotated_angles.items():
               #          f.write("{} {:10.4f} {:10.4f} 0\n".format(pos, rotated_angle, fixed_angles[pos]) )
               #     [ scale_check_vc, scale_check_hc, scale_check_M ] =   scale_check_data
               #     f.write("SC {:10.4f} {:10.4f}  {:.2f}\n".format(scale_check_vc, scale_check_hc, scale_check_M) )
               #     # hc_180deg_differences["FL"] += [ (rotated_angles["EL"], rotated_angles["WL"] - rotated_angles["EL"])   ] 
               #     # hc_180deg_differences["FR"] += [ (rotated_angles["ER"], rotated_angles["WR"] - rotated_angles["ER"])   ] 
               #     # hc_180deg_differences = 
                   
               # f.write("\n--------- hc differences --------\r")
               f.write("WIN_IDX POS-ID->ROT hc-lower(deg) hc-higher(deg) E/W-diff(arcsec)  null-LSQ-resid-std-avg\r")
               for window_index in calc_params.calculated_null_positions.keys():  # processing window indeces are in the keys
                   window = calc_params.processing_window_data[window_index]
                   rotated_angles, fixed_angles, scale_check_data, zeropos_LSQ_residstd, rotated_angles_NED = calc_params.calculated_null_positions[window_index]
                   
                   for face_dir in ["L", "R"]:
                       # if face_dir == "R": continue
                       #lower_angle_cardinal_dir  = 'E' if  rotated_angles[ 'E' + face_dir ] < rotated_angles[ 'W' + face_dir ] else 'W'
                       #higher_angle_cardinal_dir = 'W' if lower_angle_cardinal_dir == 'E' else 'E'
                       starting_cardinal_dir = 'E' if window.measurement_positions_order.index('E'+face_dir) < window.measurement_positions_order.index('W'+face_dir) else 'W'
                       ending_cardinal_dir = 'W' if starting_cardinal_dir == 'E' else 'E'
                       
                       starting_angle  = rotated_angles[ starting_cardinal_dir  + face_dir ]
                       ending_angle = rotated_angles[ ending_cardinal_dir + face_dir ]
                       
                       if ending_angle > starting_angle: 
                           passed_distance = ending_angle - starting_angle
                       else:
                           passed_distance = starting_angle - ending_angle
                       
                        
                       passed_distance_err_arcsec = (passed_distance - 180.0) * 3600.0
                       
                       if "F"+face_dir == dim_azimuth.GEODETIC_REFERENCE_FACE_DIRECTION:
                           NED_err_arcsec = (rotated_angles_NED['W'+face_dir] - rotated_angles_NED['E'+face_dir] - 180.0)  * 3600.0
                       else:
                           NED_err_arcsec = (rotated_angles_NED['E'+face_dir] - rotated_angles_NED['W'+face_dir] - 180.0)  * 3600.0
                       # resid_std_avg = np.mean( [ zeropos_LSQ_residstd[lower_angle_cardinal_dir  + face_dir] , zeropos_LSQ_residstd[lower_angle_cardinal_dir  + face_dir] ] )
                       # resid_str = "{:.3f}".format(resid_std_avg)
                       resid_str = zeropos_LSQ_residstd[starting_cardinal_dir  + face_dir] + zeropos_LSQ_residstd[ending_cardinal_dir  + face_dir]
                       f.write("{:.6f} {:} {:03} {}->{} HOR-ROT:  {:8.4f} -> {:8.4f} mismatch: {:+10.4f} arcsec |  NED-mismatch: {:+10.4f} arcsec | resid:{}\n".format(window.first_M_timelabel, mu.exceltime_2_datetimestr(window.first_M_timelabel), window_index, starting_cardinal_dir+face_dir, ending_cardinal_dir+face_dir, starting_angle, ending_angle, passed_distance_err_arcsec, NED_err_arcsec, resid_str ) )
                
               saved_files += [calc_params.calculated_null_positions_file]
               
    if nr_of_all_caclvars_valid_results  ==  0:
         print("There is no valid result at all!")
         warn_logs += DIM_data_loader.get_warnings()
         return None, None, saved_files, warn_logs

    if enable_post_plots:

        for icalcvar, calc_params in enumerate(calc_params_list):
            
            if calc_params.nr_of_valid_results == 0: continue
            # Creating a new plot for pos_id slope
            pos_id_list = np.unique(calc_params.ra_measpoints["pos_ID"])
            pos_id_avg_slope_orig = np.empty( len(pos_id_list) )
            pos_id_avg_slope_corrected_by_variometer = np.empty( len(pos_id_list) )

            for i,pos_id in enumerate(pos_id_list):
                pos_id_avg_slope_orig[i] = np.mean(  calc_params.ra_measpoints[ calc_params.ra_measpoints["pos_ID"]==pos_id    ]["M_meas_slope"]   )
                # pos_id_avg_slope_corrected_by_variometer[i] = pos_id_avg_slope_orig[i] - np.mean( calc_params.ra_measpoints[ calc_params.ra_measpoints["pos_ID"]==pos_id    ]["M_slope_from_variometer"] )

            plot_utils.plot_xy( pos_id_list, [pos_id_avg_slope_orig, pos_id_avg_slope_corrected_by_variometer], "pos-id", "Slope (nT/s)", ["x", [".","-"][1] ] , color = ["blue", "red"], label = ["oig", "corrected"], close_plot = False )
            for i,pos_id in enumerate(pos_id_list):
                plt.gca().annotate( "{:}".format(pos_id)  , (pos_id, pos_id_avg_slope_corrected_by_variometer[i]))
            plt.show()


            # Calculating and plotting residual averages for each pos-id.
            # The average residual can depend, for example, on the eccentricity error
            if len(calc_params.resid_per_pos_id_avgsum)>0:
                resid_per_pos_id_avg = {}
                for pos_id in calc_params.resid_per_pos_id_avgsum:
                    resid_per_pos_id_avg[pos_id] = calc_params.resid_per_pos_id_avgsum[pos_id] / calc_params.resid_per_pos_id_count[pos_id]
                pos_id_resid_array = np.array(list(resid_per_pos_id_avg.items())) # items includes keys and values
                plot_utils.plot_xy(pos_id_resid_array[:,0],pos_id_resid_array[:,1],"pos-ID of calc-var-{}".format(icalcvar), "avg-residual [nT]", close_plot = False)
                for i in range( len(pos_id_resid_array)) :
                    plt.gca().annotate( "{:}".format(int(pos_id_resid_array[i,0]))  , (pos_id_resid_array[i,0], pos_id_resid_array[i,1]))
                plt.show()

    if enable_post_plots and (nr_of_all_caclvars_valid_results>0):

        x_axis_compname = "datetime-mid"  # "cycle_nr" or "datetime"

        # Reading in reference data (1-min variation) for each calc variation
        # dim_measurement_loader already read and stored the reference data at start, anyway, we read it out once again,
        # as here we need some more data at the beginning and end

        calc_params_list_minute_ref_data_timelabels = []
        calc_params_list_minute_ref_data_comps = [] #

        if (mag_obs_dataset is not None) and mag_obs_dataset.has_variometer_data() and mag_obs_dataset.has_F_data():

            for calc_params_index, calc_params in enumerate( calc_params_list ) :
                
                if calc_params_list[calc_params_index].nr_of_valid_results == 0:
                    calc_params_list_minute_ref_data_timelabels += [ [] ]
                    calc_params_list_minute_ref_data_comps      += [ [] ]
                    continue

                # First time label of first processing window
                t0 = calc_params.result_array[0]["datetime-start"]
                # Finding the last valid time label of the last processing window:
                t1 = calc_params.result_array["datetime-end"]
                # Extracting filled timelabels sub-array (timelabels not having nan values)
                t1 = t1[ np.invert( np.isnan(t1) ) ]
                t1 = t1[-1]

                assert isinstance(t0,float) and isinstance(t1,float)

                ref_timelabel_start = t0  - 0*(t1-t0)/20  # Starting earlier than first time label
                ref_timelabel_end   = t1  + 0*(t1-t0)/20  # Ending later than last time label
                ref_minute_timelabels = dat_file_reader.get_minute_timelabels(ref_timelabel_start, ref_timelabel_end)

                ref_data_minute_variation, missing_record_indeces = mag_obs_dataset.load_data(True, False, ref_minute_timelabels,  result_data_type = dat_file_reader.file_read_result_data_type._1min_arithm_avg, in_case_of_missing_error = "none" ) # data_descr = "1mindata-for-plot", 
                
                # Not used parameter: measurement_sites.F_default

                calc_params_list_minute_ref_data_timelabels += [ ref_minute_timelabels ]
                calc_params_list_minute_ref_data_comps      += [ ref_data_minute_variation ]


        # concatenated_plot_y_components = concatenated_plot_y_components + ["ae1_hc_ampl", "ae1_hc_phase","ae1_vc_ampl", "ae1_vc_phase", ... ]

        for plotted_comp_names in concatenated_plot_y_components:

            # We plot one or more component in a plot - we convert simple comp name definitions into a one-element list
            plotted_comp_names = mu.ensure_list(plotted_comp_names)

            # Special plot:
            # Different symbols for different calc. variations ('plot_symbol').
            # Different colors for
            #     - different methods ('plot_color').
            #     - If we plot more components in one plot, then for different components (deltaX/deltaY/deltaZ)

            nof_components_plotted_together = len(plotted_comp_names)
            y_axis_label = ""

            # 'plotted_comp_names' can have more component together ('nof_components_plotted_together'>1)
            # Now we plot each:
            for plotted_comp_index,plotted_comp_name in enumerate(plotted_comp_names):

                print("Plotting ""{}""...".format(plotted_comp_name))
                comp_unit = None
                for calc_params_index, calc_params in enumerate( calc_params_list ) :
                    
                    if calc_params.nr_of_valid_results == 0: continue

                    if not plotted_comp_name in calc_params.results_array.dtype.names:
                        # Not all component exist in all calc. variations
                        continue

                    if comp_unit is None: comp_unit = result_array_comps_params[plotted_comp_name]["unit"]
                    plot_symbol = calc_params_plot_symbols[calc_params_index]

                    if glob_test_mode in [ "temp-dependence", "tempchange-dependence" ]:
                        plot_symbol='--bo'

                    # Getting results for the plotted component of 'plotted_comp_name' and for the actual calc. var:

                    ra_result = calc_params.results_array[plotted_comp_name]
                    ra_x_comp = calc_params.results_array[x_axis_compname]

                    # Plotting points from different methods in different colots:
                    ra_main_program_type_idx = calc_params.results_array["main_program_type_idx"]

                    for main_program_type_idx, program_main_type_str in enumerate(DIM_data_loader.program_maintype_str_list):
                         x = ra_x_comp[  ra_main_program_type_idx ==  main_program_type_idx ]
                         y = ra_result[  ra_main_program_type_idx ==  main_program_type_idx ]
                         if (len(x)>0):
                             plot_color = plot_utils.rainbow_color_get( plotted_comp_index *  len(DIM_data_loader.program_maintype_str_list) +  main_program_type_idx, len(DIM_data_loader.program_maintype_str_list) * nof_components_plotted_together )
                             plot_label = calc_params.name
                             if len(DIM_data_loader.program_maintype_str_list) > 1: plot_label +=" / "+meas_method_str # .decode('UTF-8')
                             if nof_components_plotted_together>1: plot_label=plotted_comp_name+' '+plot_label
                             # print("  with plot label of ""{}""".format(plot_label))
                             if x_axis_compname == "datetime-mid":
                                 plot_utils.plot_datetime_y(x, [y], plotted_comp_name, label=[plot_label], color = [plot_color], b_close_plot = False, fmt=plot_symbol)
                             else:
                                 plt.plot( x, y, plot_symbol, color = plot_color, label = plot_label  )

                    # Adding reference plots to the already existing plot:
                    if len(calc_params_list_minute_ref_data_comps)>0 and  ( plotted_comp_name in calc_params_list_minute_ref_data_comps[calc_params_index].dtype.names ):
                        # plotting DIDD'S I/D line to the determined I/D
                        assert x_axis_compname == "datetime-mid"
                        plot_utils.plot_datetime_y( calc_params_list_minute_ref_data_timelabels[calc_params_index], [ calc_params_list_minute_ref_data_comps[calc_params_index][plotted_comp_name] ], plotted_comp_name, fmt='-', color=['red'], label=['Ref-'+plotted_comp_name], b_close_plot = False)

                if plotted_comp_index>0: y_axis_label +=  " "
                assert not comp_unit is None, "No unit for plotted comp \"{}\"".format(plotted_comp_name)
                y_axis_label += plotted_comp_name + " [" + comp_unit + "]"

            plt.xlabel(x_axis_compname)
            plt.ylabel(y_axis_label)
            if "X0" in plotted_comp_names:
                plt.gca().set_yticks( [-5, -2.5, 0, +2.5, +5] )
                plt.gca().grid(axis='y', linestyle="--", linewidth=0.5, color = 'red' )

            if b_more_calc_params or ( len(DIM_data_loader.program_maintype_str_list) > 1) or (nof_components_plotted_together>1):
                # plt.legend(framealpha=0.5, frameon=True, bbox_to_anchor=(1.1, 1.05));
                # plt.legend(framealpha=0.5, frameon=True, loc = 'best');
                pass

            if glob_test_mode in [ "temp-dependence", "tempchange-dependence" ]:
                if (x_axis_compname == "datetime-mid") and (plotted_comp_name in ['I-diff','D-diff','delta','epsilon','M_scale','M_offset','Z-diff']):
                    temp_axes = plt.gca().twinx()
                    if glob_test_mode == "temp-dependence":
                        temp_axes.set_ylabel("Temp (°C)")
                        assert False
                        temp_axes.plot(  mu.exceltime_to_matplotlib_time(ra_x_comp),  calc_params[0].result_array.comp_array_get("Temp"), color="tab:brown", linewidth=4.0, linestyle='dashed' )
                    elif glob_test_mode == "tempchange-dependence":
                        temp_axes.set_ylabel("Temp-change (°C)")
                        temp_axes.plot(  mu.exceltime_to_matplotlib_time(ra_x_comp),  calc_params[0].result_array.comp_array_get("Temp-change"), color="tab:brown", linewidth=4.0, linestyle='dashed' )
                    else:
                        assert False


            plt.show()

    # T = result_array[0].comp_array_get("Temp")
    # def temperature_dependence_plot(comp_name):
    #     plt.plot(T,  result_array[0].comp_array_get(comp_name), "*")
    #     conv_A, conv_M = mu.LSQ( T ,  result_array[0].comp_array_get(comp_name) )
    #     T_linspace= np.linspace(min(T), max(T), 100 )
    #     Y = conv_A + conv_M * T_linspace
    #     plt.plot(T_linspace,Y, "-", color = "red")
    #     plt.xlabel('Temperature (°C)')
    #     plt.ylabel(comp_name+" ["+result_array[0].comp_unit_get(comp_name)+"]")
    #     plt.title(comp_name+ " vs. temperature");
    #     plt.show()
    #     print( f"M_scale: {conv_A=}  {conv_M=}" )

    # temperature_dependence_plot("M_scale")
    # temperature_dependence_plot("M_offset")
    # temperature_dependence_plot("epsilon")
    # temperature_dependence_plot("delta")

    if do_tick_mod_4_analysis:
        # ADS digital angle reader introduced error estimation (deprecated)
        tick_mod_4_symbols = [ 'o',  'x', '.', 'v' ]
        for icalcvar in range(nof_calc_params):
            print("--- tick-mod-4 anaylsis of paramss {} ---".format(icalcvar))
            if len(ra_measpoints[icalcvar])==0:
                cu.warning("No result for calc-var #{}".format(icalcvar))
                continue
            tick_mod_4_bins_hc = [ [],[],[],[] ] 
            tick_mod_4_bins_vc = [ [],[],[],[] ] 
            
            for r in ra_measpoints[icalcvar]:
                index = r['hc_tick_mod_4']
                tick_mod_4_bins_hc[ index ] += [ r['M_mean_removed_resid_nT'] ]
                index = r['vc_tick_mod_4']
                tick_mod_4_bins_vc[ index ] += [ r['M_mean_removed_resid_nT'] ]
            
            print("-- hc --")
            for index in [0,1,2,3]:
               print("hc-tick-mod-{:}: count:{:}({:.1%})  avg:{:.2f}".format(index, len(tick_mod_4_bins_hc[index]), len(tick_mod_4_bins_hc[index]) / len(ra_measpoints[icalcvar]), np.mean( tick_mod_4_bins_hc[index] )   ) )
            print("-- vc --")
            for index in [0,1,2,3]:
               print("vc-tick-mod-{:}: count:{:}({:.1%})  avg:{:.2f}".format(index, len(tick_mod_4_bins_vc[index]), len(tick_mod_4_bins_vc[index]) / len(ra_measpoints[icalcvar]), np.mean( tick_mod_4_bins_vc[index] )   ) )
                 
    #             x = []
    #             y = []
    #             for sample_nr in range(len(hc_tick_raw)):
    #                 if ( round(hc_tick_raw[sample_nr]) % 4) == tick_mod_4:
    #                     x = x + [sample_nr]
    #                     y = y + [ best_res_residuals[sample_nr] ]
    #         if len(x) > 0:
    #             plt.plot(x, y, tick_mod_4_symbols[tick_mod_4], label = "tick % 4 = {}".format(tick_mod_4))
    #             plt.xlabel('sample number')
    #             plt.ylabel('residual of the best solution (nT)')
    #             plt.legend(framealpha=1, frameon=True)
    #             plt.show()
        
    
    # if len( ra_measpoints[0] ) >4 :
    #   for icalcvar,calcvar in enumerate(calc_params_list):
    #       optimize_eccentricity_errors_based_on_resid_array(icalcvar, calcvar, ra_measpoints[icalcvar] )
    
    warn_logs += DIM_data_loader.get_warnings()
    return (min_result_time, max_result_time), None, saved_files, warn_logs


# Test function if time shift can be estimated from residuals
def dim_time_shift_estimation( dat_files_path, filenames, DIM_instr, calc_params_list, mag_obs_dataset ):

    assert len(calc_params_list)==1
    enable_cycle_plots = False
    enable_cycle_tables = False
    enable_post_plots = False
    enable_post_process = False 
    b_remove_means = False 
    resid_avg_list = []
    tested_time_shifts_sec = np.arange(-0.5,  +2.0, 0.25) 
    tested_time_shifts_sec = np.arange(-0.7,  0.0,  0.02) 
    for i,time_shift_sec in enumerate(tested_time_shifts_sec):
        
        variometer_time_shift_sec = -0.39
        F_time_shift_sec = time_shift_sec
        
        ra_cycles, ra_measpoints = dim_measurement_files_process(dat_files_path, filenames, DIM_instr, calc_params_list, enable_cycle_plots, enable_cycle_tables, enable_post_plots, enable_post_process, mag_obs_dataset = mag_obs_dataset, b_remove_means = b_remove_means, variometer_time_shift_sec = variometer_time_shift_sec, F_time_shift_sec = F_time_shift_sec, processing_mode = "PER-CYCLES", meas_pos = None )
        # resid_std = np.nanmean(ra_cycles[0]["resid_std"])
        resid_std = np.std( ra_measpoints[0]["M_mean_removed_resid_nT"] )  
        # vagy:  np.std( ra_measpoints[0]["M_raw_resid_nT"] )
        print("========================== {}/{} avg resid-stddev: {} nT =========================".format(i+1,len(tested_time_shifts_sec), resid_std))
        resid_avg_list += [ resid_std ] 
        
    plot_utils.plot_xy(tested_time_shifts_sec, resid_avg_list, "time-shift (sec)", "resid-avg (nT)") 
