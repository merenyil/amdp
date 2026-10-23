# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 20:44:03 2025

@author: László Merényi 

Variables / functions for general geomagnetic field calculations
and for management of geomagnetic observatory specific data.

For DIM/ADS specific staff, see dim_utils.py and absolute_magnetic_data_processing.py files
"""

import numpy as np
import math
import yaml
from geomag_utils_preferences import *
import misc_utils as mu
import amdp_base_formulas


BASELINE_CALCULATION_APPROACHES = ("auto", "none", "D0-angle", "Y0-nT", "dIdD"  )
COMPNAME_VARIOMETER_H = 'H'
COMPNAME_VARIOMETER_D = 'D'
COMPNAME_VARIOMETER_Z = 'Z'
COMPNAME_TOTAL_H      = 'h'
COMPNAME_DECL_ANGLE   = 'd'
COMPNAME_INCL_ANGLE   = 'i'
COMPNAME_DECL_INTENSITY='D'
COMPNAME_INCL_INTENSITY='I'

COMPNAME__BASELINE__HX_SENSOR = "H0"
COMPNAME__BASELINE__DY_SENSOR = "D0_nT"
COMPNAME__BASELINE__D_ANGLE   = "D0_rad"  # 'decbas'
COMPNAME__BASELINE__Z_SENSOR  = "Z0"
# COMPNAME__BASELINE__DECBAS   = "decas"  # Sensor orientation in radian
COMPNAME__F_PILLAR_DIFF = "F_pillar_diff"


IAGA2002_DATA_TYPE__AUTOMATIC    = "A"
IAGA2002_DATA_TYPE__DEFINITIVE   = "d"
IAGA2002_DATA_TYPE__PROVISIONAL  = "p"
IAGA2002_DATA_TYPE__QUASI_DEFINITIVE = "q"
IAGA2002_DATA_TYPE__VARIATION    = "v"

PILLAR_NAME__ALL_PILLARS = "*"

MASK_PIER1 = "${PIER}"   
MASK_PIER2 = "${PILLAR}"   
MASK_SENSORPOS = "${SP}"
MASK_THEO_FACE_POSITION = "${FP}"
MASK_REFPOS    = "${RP}"

measurement_site_list = {}

# # Plot colors and label names if we want to plot XYZF data
# XYZF_plot_colors = ("black","blue","green","red")
# XYZF_names  = ('X','Y','Z','F')


def threecompmag_to_scalarfield_model(model_type, X,Y,Z, model_params):
    assert len(model_params) == 9
    X_offset = model_params[0]
    Y_offset = model_params[1]
    Z_offset = model_params[2]
    X_scale = model_params[3]
    Y_scale = model_params[4]
    Z_scale = model_params[5]
    alpha = model_params[6]
    beta = model_params[7]
    gamma = model_params[8]


    if model_type=="ML-model":
        m_X_corr = (X - X_offset)*X_scale
        m_Y_corr = (Y - Y_offset)*Y_scale
        m_Z_corr = (Z - Z_offset)*Z_scale
        # http://xraylab.chem.purdue.edu/lesson2.ppt
        F_corr =  np.sqrt(np.power(m_X_corr,2) + np.power(m_Y_corr,2) + np.power(m_Z_corr,2) + 2*m_X_corr*m_Y_corr*np.cos(np.pi/2+gamma) + 2*m_X_corr*m_Z_corr*np.cos(np.pi/2+beta) + 2*m_Y_corr*m_Z_corr*np.cos(np.pi/2+alpha) )
    elif model_type=="ADAM-model":
        # Ádám modell:
        m_x_corr = X_scale*(X - X_offset) + gamma*(Y - Y_offset) + beta*(Z - Z_offset)
        m_y_corr = Y_scale*(Y - Y_offset) + alpha*(Z - Z_offset)
        m_z_corr = Z_scale*(Z - Z_offset)
        F_corr = np.sqrt( np.square(m_x_corr) + np.square(m_y_corr) + np.square(m_z_corr) )
    else:
        assert False

    return F_corr

# Converts exactly XY oriented variometer data into Declination angle(s).
# Orientation error and orthogonality errors are not taken into account.
def XY_2_Decl(X,Y):
    return np.arctan2(Y,X)
    # Other possible equation: cos(D)=X/H

# Converts HDZ oriented variometer recording into Declination angle.
# Orthogonality errors are not taken into account.
def NE_2_Decl(N,E,decbas):  # ~HD_2_Decl~
    return np.arctan( E / N) + decbas  # Or: np.arctan2( E, N )

# 'H' is full horizontal component, 'E' is the 'D' axis measurement
def HE_2_Decl(H,E,decbas):  #
    return np.arcsin( E / H) + decbas

def totalH_Z__2__Incl(totalH,Z):
    F = np.sqrt( np.square(totalH) + np.square(Z))  # Rossz volt a számítás 2026.07.29-ig, mert a totalH nem volt négyzetre emelve!!!
    # print("---------------")
    # print(F)

    r = totalH / F
    assert np.all(r >= -1)
    assert np.all(r <= +1), "{}".format(r)

    I = np.arccos( r )
    return I

def XYZ_2_Incl(X,Y,Z):
    assert False, "Consider using 'totalH_Z__2__Incl"
    XX_YY = X*X + Y*Y
    H = np.sqrt(XX_YY)
    F = np.sqrt(XX_YY + Z*Z)
    I1 = np.arccos( H / F )

    # Other possible expression: sin(I) = Z/F (?)
    # I2 = np.arcsin( Z / F)
    # This check fails if data set has NAN's
    # if isinstance(I2,float):
    #     assert abs(I1-I2) < 1e-3
    # else:
    #    assert np.max(np.abs(I1-I2)) < 1e-3, f"{X=} {Y=} {Z=}"
    return I1

def IncleDecl_2_XYZ(Incl,Decl,F):
    Z = F * np.sin(Incl)
    H = F * np.cos(Incl)
    X = np.cos(Decl) * H
    Y = np.sin(Decl) * H
    return X,Y,Z

def XY_2_H(X,Y):
    return np.sqrt( np.square(X) + np.square(Y) )


def XYZ_2_F(X,Y,Z):
    return np.sqrt( np.square(X) + np.square(Y) + np.square(Z) )

# 'H_or_X' and 'D_or_Y' and 'Z' must be baseline corrected axis components in nT units
def calc_D_I_angles(decbas, H_or_X, D_or_Y, Z):

   D_angle = NE_2_Decl(H_or_X, D_or_Y,decbas)

   # if sensor_orientation=="HDZ":
   #    variometer_orientation_to_measurement_time_declination = np.arctan2(D_or_Y, H_or_X)
   #    D_angle = decbas + variometer_orientation_to_measurement_time_declination
   #elif sensor_orientation=="XYZ":
   #    # We have to return the nT values that ADS should measure in hc_theo/vc_theo positions, and also fits to F
   #    	variometer_orientation_to_measurement_time_declination = np.arctan2( D_or_Y , H_or_X ) # Or: np.arctan( Y / Y )
   #        D_angle = decbas + variometer_orientation_to_measurement_time_declination

   Horiz = np.sqrt( np.square(H_or_X) +  np.square(D_or_Y) )
   I_angle = np.arctan2( Z, Horiz )



   return D_angle,I_angle

def HD_2_XY(H,D,decbas):
    X = H * np.cos(decbas) - D * np.sin(decbas)
    Y = H * np.sin(decbas) + D * np.cos(decbas)
    return X,Y

# Converts variometer data (XYZ or HDZ) into (preliminary) absolute values at absolute pillar.
# 'raw_varimeter_data' must be a structured array
# 'baseline_data' is a dictionary, supplying conversion rules (baselines)
def variometer_to_abs_data(calc_approach, mag_obs_data_def, calc_params_common, raw_varimeter_data, baseline_data, allow_empty_baseline = True ):

    # result = np.empty( (len(raw_varimeter_data)), dtype=[ (COMPNAME_INCL_ANGLE,float),  (COMPNAME_DECL_ANGLE,float), ('F',float) ] )

    # F_pillar_diff = pillar_diff_def.F_pillar_diff_variometer

    assert calc_approach in BASELINE_CALCULATION_APPROACHES[1:]


    if calc_approach == "Y0-nT": # mag_obs_data_def.instr_frame()[:4] == "XYZ":

        XH_compname = mag_obs_data_def.usage_frame[0] # raw_varimeter_data.dtype.names[0]
        assert XH_compname in ['X','H']
        YD_compname = mag_obs_data_def.usage_frame[1] # raw_varimeter_data.dtype.names[1]
        assert YD_compname in ['Y','D']

        if (baseline_data is not None) and (COMPNAME__BASELINE__HX_SENSOR in baseline_data.keys()):
            XH_basline = baseline_data[COMPNAME__BASELINE__HX_SENSOR]
        elif allow_empty_baseline:
            XH_basline  = 0.0
        else:
            assert False, "{}".format(baseline_data)

        if (baseline_data is not None) and (COMPNAME__BASELINE__DY_SENSOR in baseline_data.keys()):
            YD_baseline  = baseline_data[COMPNAME__BASELINE__DY_SENSOR]
        elif allow_empty_baseline:
            YD_baseline  = 0.0
        else:
            assert False, "{}".format(baseline_data)

        if (baseline_data is not None) and (COMPNAME__BASELINE__Z_SENSOR in baseline_data.keys()):
            Z_baseline  = baseline_data[COMPNAME__BASELINE__Z_SENSOR]
        elif allow_empty_baseline:
            Z_baseline  = 0.0
        else:
            assert False, "{}".format(baseline_data)


        XH_abs = raw_varimeter_data[XH_compname] + XH_basline # X0 or H0
        YD_abs = raw_varimeter_data[YD_compname] + YD_baseline # Y0 or D0_nT
        H_abs = np.sqrt( np.square(XH_abs) + np.square(YD_abs) )
        Z_abs = raw_varimeter_data['Z'] + Z_baseline
        F_abs = np.sqrt( np.square(XH_abs) + np.square(YD_abs) + np.square(Z_abs) ) + calc_params_common.F_pillar_diff_variometer
        I_abs = np.arctan2(Z_abs, H_abs)
        # I_abs_check = np.arcsin(Z_abs/F_abs)
        # I_diff = I_abs_check-I_abs
        # I_diff= I_diff[~np.isnan(I_diff)]
        # assert np.max(np.abs(I_diff)) < np.radians( 10/3600. ), "{:.1f} arcsec".format(np.degrees(np.max(np.abs(I_abs_check-I_abs)))*3600.0)
        D_abs = np.arctan2(YD_abs, XH_abs) + mag_obs_data_def.get_decbas()  # ToDo: check if decbas must be added or removed
        
    elif calc_approach == "D0-angle": #  mag_obs_data_def.instr_frame()[:4]=="HDZ":

        H_compname = mag_obs_data_def.usage_frame[0]  # raw_varimeter_data.dtype.names[0]
        D_compname = mag_obs_data_def.usage_frame[1] # raw_varimeter_data.dtype.names[1]
        assert H_compname in ['H','X']  # COMPNAME_VARIOMETER_H
        assert D_compname in ['D','Y']  # COMPNAME_VARIOMETER_D
        H_var = raw_varimeter_data[H_compname]
        E_var = raw_varimeter_data[D_compname]
        Z_var = raw_varimeter_data['Z']
        decbas = None
        D_baseline = None
        baseline_system = None

        if (baseline_data is None) and allow_empty_baseline:
            baseline_data = {}
            baseline_data[COMPNAME__BASELINE__HX_SENSOR] = 0.0
            baseline_data[COMPNAME__BASELINE__D_ANGLE] = 0.0
            baseline_data[COMPNAME__BASELINE__Z_SENSOR]  = 0.0

        for key,value in baseline_data.items():
            if   key==COMPNAME__BASELINE__HX_SENSOR:
                H_abs = H_var + value
                H_var = None
            elif key==COMPNAME__BASELINE__DY_SENSOR:
                assert baseline_system is None
                baseline_system = "AXES"
                E_abs = E_var + value
                E_var = None
            elif key==COMPNAME__BASELINE__Z_SENSOR:
                Z_abs = Z_var + value
                Z_var = None
            #elif key==COMPNAME__BASELINE__DECBAS:
            #    assert decbas is None # 'decbas' should be set only once
            #    decbas = value
            elif key == COMPNAME__F_PILLAR_DIFF:
                assert False, "Use other pillar_diff_def parameter instead"
                F_pillar_diff = value
            elif key == COMPNAME__BASELINE__D_ANGLE:
                assert baseline_system is None
                baseline_system = "DECL_ANGLE"
                assert D_baseline is None
                D_baseline = value
                E_abs = E_var
                E_var = None
            elif key in [ "mag_sensor", "I", "D" ]:
                pass # 'baseline_data' can have this field, if it is a 'last_result' of a calculation variation
            else:
                assert False, "Unknown baseline parameter: \"{}\"".format(key)
        if baseline_system == "AXES":
            assert decbas is not None
            assert D_baseline is None
            D_abs,I_abs = calc_D_I_angles(decbas, H_abs, E_abs, Z_abs)
        else:
            assert decbas is None
            if (D_baseline is None): return None, None, None
            D_abs, I_abs = calc_D_I_angles(D_baseline, H_abs, E_abs, Z_abs)

        # QQQ: pillar difference alkalmazás sorrendje!
        I_abs += calc_params_common.dI_pillar_diff
        D_abs += calc_params_common.dD_pillar_diff
        F_abs = np.sqrt( np.square(H_abs) + np.square(E_abs) + np.square(Z_abs) ) + calc_params_common.F_pillar_diff_variometer

    else:
        assert False, calc_approach
        
    result = magfield_class( I = I_abs, D = D_abs, F = F_abs)
    return result # F_abs,I_abs,D_abs        
        
# A magnetic field vector model, based on known I/D and F values.
class magfield_class:

    def __init__(self, meas_pos = None, I=None, D=None, F=None ):
        self.reference_coord_system = 'theo'  # X/Y/Z are valid for this system
        if meas_pos is None :
            # assert default_meas_pos is not None # The default measurement position must be selected prior calling this
            # meas_pos = default_meas_pos
            self.F = F
            self.D = D
            self.I = I
        else:
            self.F = F if F is not None else meas_pos.F_typical
            self.D = D if D is not None else meas_pos.D_typical
            self.I = I if I is not None else meas_pos.I_typical
        self.last_set_X = np.nan
        self.last_set_Y = np.nan
        self.last_set_Z = np.nan


    def X(self):
        return self.H() * np.cos(self.D)
        # _X,_Y,_Z = ID_2_XYZ(self.I,self.D,self.F)
        # return _X

    def Y(self):
        return self.H() * np.sin(self.D)
        # _X,_Y,_Z = ID_2_XYZ(self.I,self.D,self.F)
        # return _Y

    def H(self):
        return self.F * np.cos(self.I)

    def calc_horiz_comp(self, NED_direction):
        return np.cos( self.D - NED_direction ) * self.H()

    def Z(self):
        # return self.F * np.cos(self.I - np.pi/2)
        return self.F * np.sin(self.I)

    def set_from_XYZ(self,X,Y,Z):
        self.last_set_X = X
        self.last_set_Y = Y
        self.last_set_Z = Z
        self.D = XY_2_Decl(X,Y)
        self.I = totalH_Z__2__Incl(XY_2_H(X,Y), Z)
        self.F = math.sqrt(X*X + Y*Y + Z*Z)


def get_decimal_digit_for_unit(comp_unit):
    # Auto determining number of decimal digits
    if comp_unit=="": return 0
    if comp_unit=="nT": return nT_PREC
    if comp_unit in  ['°','deg','gon']:  return INCL_DECL_DEG_PREC
    if comp_unit in  ['min']:      return 1
    if comp_unit in  ['"','arcsec', "'", 'arcmin']: return 1
    if comp_unit=="rad": return 6
    if comp_unit=="tick": return 1
    if comp_unit=="K": return 2 # Kelvin (as temperature)
    assert False, "Not supported unit: {}".format(comp_unit)

        

class measurement_pos_class:

    def __init__(self, station_code, pillar_def = "*", set_as_default = True):

        self.station_code = station_code
        if station_code in measurement_site_list.keys():
            sitedef = measurement_site_list[station_code]
            self.I_typical  = math.radians(sitedef['I_deg'])
            self.D_typical  = math.radians(sitedef['D_deg'])
            self.F_typical  = sitedef['F']
            self.sites_pillars    = sitedef['piers']
            self.azimuth_list   = sitedef['azimuth-list']
        else:
            assert False, f'Station "{station_code}" not found in stations list!'
            self.I_typical = math.radians(45.0)
            self.D_typical = 0.0
            self.F_typical = 50000.0
            self.sites_pillars = mu.ensure_list(pillar_def)

        self.pillars = mu.ensure_list(pillar_def)
        
        if pillar_def == PILLAR_NAME__ALL_PILLARS:
            pass
        else:
            for pillar_check in mu.ensure_list(self.pillars):
                assert pillar_check in self.sites_pillars, "Pier \"{}\" is not found in the sites/observatory's piers list".format(pillar_check)

        self.D_min = self.D_typical - math.radians(4)
        self.D_max = self.D_typical + math.radians(4)
        self.I_min = self.I_typical - math.radians(4)
        self.I_max = self.I_typical + math.radians(4)


        #if set_as_default:
        #    global default_meas_pos
        #    default_meas_pos = self

    #def refpos_const_data_source_records(self):

    #    return measurement_site_list[self.station_code]["azi-data"]




# Builds 'measurement_site_list' and 'refpos_const_data_list'
def pillar_data_init(_measurement_sites_yaml_path):
    
    with open(_measurement_sites_yaml_path, "r", encoding="utf-8") as _f:
        _measurement_sites_yaml = yaml.safe_load(_f)
        
    for _site_code, _site_def in _measurement_sites_yaml.items():
        
        assert not _site_code is measurement_site_list.keys()
        
        _site_def = dict(_site_def)  # shallow copy, so popping "azimuth-marks" doesn't touch the loaded YAML data
        if "azimuth-marks" in _site_def:
            azi_list = _site_def.pop("azimuth-marks")
        else:
            azi_list = []
        measurement_site_list[_site_code] = _site_def                            
        
        # Building list of azimuth marks
        # assert not _site_code in refpos_const_data_list.keys() 
        measurement_site_list[_site_code]["azimuth-list"] = []
        for azi_def in azi_list:
            azi_def =  azi_def.replace(',', ' ')            
            azi_def =  azi_def.split()
            azi_dict = {}            
            azi_dict["pier"] = azi_def[0]
            azi_dict["azimark-name"] = azi_def[1]
            azi_dict["hc-dir"] = math.radians(float(azi_def[2]))
            if len(azi_def)>3:
                azi_dict["vc-dir"] = math.radians(float(azi_def[3]))
            else:
                azi_dict["vc-dir"] = None
            
            measurement_site_list[_site_code]["azimuth-list"] += [azi_dict]
            
                 
                 
# filename can be a single filename or a list/tuple of filenames                 
def replace_masks_in_filename( filename, meas_pos, pillar_index = 0 ):
    if isinstance(filename, (list,tuple)) :
        for i in range( len(filename)  ) :
            filename[i] = replace_masks_in_filename( filename[i], meas_pos )
    else:
        assert isinstance(meas_pos, measurement_pos_class )
        if MASK_PIER1 in filename:
            assert len(meas_pos.pillars) == 1
            filename = filename.replace(MASK_PIER1, meas_pos.pillars[pillar_index] ) 
        if MASK_PIER2 in filename:
            assert len(meas_pos.pillars) == 1
            filename = filename.replace(MASK_PIER2, meas_pos.pillars[pillar_index] ) 
        
    return filename
        

        
def M_recordings_time_reduction( reduct_to_index, M, hc_NED, vc_NED, abs_magfield, sensor, correction = 1.0 ):
    assert reduct_to_index < len(M)

    M_variometer_at_ref_time_field = amdp_base_formulas.dim_instrument_model(abs_magfield.F[reduct_to_index], vc_NED, hc_NED, abs_magfield.I[reduct_to_index], abs_magfield.D[reduct_to_index], sensor.epsilon(), sensor.delta(), sensor.offset, sensor.scale)
    M_variometer_at_changing_field = amdp_base_formulas.dim_instrument_model(abs_magfield.F,                  vc_NED, hc_NED, abs_magfield.I,                  abs_magfield.D,                  sensor.epsilon(), sensor.delta(), sensor.offset, sensor.scale)
    M_corr = M_variometer_at_changing_field - M_variometer_at_ref_time_field
    
    return M - M_corr * correction
    


