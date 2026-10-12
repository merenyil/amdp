"""
dim_utils.py

Variables and functions for DIM (including ADS) measurement specific calcualtions,
used for baseline calculations. Baseline value calculation are not here.

See also 'absolute_magnetic_data_processing.py'
@author: László Merényi

"""

import math
import dim_instr
import amdp_base_formulas
import dim_azimuth
import numpy as np
import geomag_utils as gu
import angle_utils
from enum import Enum


class telescope_face_pos_enum(Enum):
    FL = 'L'   # Or: "direct position"
    FR = 'R'   # Or: "reverse position"
    unknown = '?'


UNKNOWN_CARDINAL_DIRECTION = '?'
TELESCOPE_CARDINAL_DIRECTIONS = "NESW"  # North/East/South/West

SCALECHECK_POS_ID = "SC"
SCALECHECK_POS_NAME = "Scale_check"
# std_NM_meas_part_names   = ( "ED", "WD", "EU", "WU", "NU", "SD", "ND", "SU", SCALECHECK_POS_NAME)
abs_meas_positions_ids     = ( "EL", "WL", "ER", "WR", "NL", "SR", "NR", "SL", SCALECHECK_POS_ID)
abs_meas_positions_names   = ( "E/FL", "W/FL", "E/FR", "W/FR", "N/FL", "S/FR", "N/FR", "S/FL", SCALECHECK_POS_NAME)
abs_meas_positions_declincl_indeces = (0,1,2,3,4,5,6,7)  # slice(0,8)
abs_meas_positions_decl_indeces = (0,1,2,3)  #  slice(0,4)
abs_meas_positions_incl_indeces = (4,5,6,7)  # slice(4,8)
abs_meas_positions_declincl_slice = slice(0,8)
abs_meas_positions_decl_slice =  slice(0,4)
abs_meas_positions_incl_slice =  slice(4,8)


conventional_DIM_geodetic_pos_2_sensor_pos = { 'L': 'U', 'R': 'D' }  # Face left is "Geodetic position I", and sensor is UP in this position


# F_DATA_FLAGS = [gu.F_DATA_SOURCES.VARIOMETER_DATA, gu.F_DATA_SOURCES.SITES_DEFAULT_F, gu.F_DATA_SOURCES.DIM_DATA]
# MAGDATA_DEF__PILLAR_DIF_INDEX = 1   # ref_data_def  and varcorr_data_def is a list/tuple, in where 'MAGDATA_DEF__PILLAR_DIF_INDEX' tells which element is the pillar difference


# This function tries to split the ADS's measurement points into main measurement parts.
# It is used for residual plots, where we marker of the different parts of nullmethod are different.
# Deprecated from 2025-05. Use 'separate_DIM_measurement_positions' instead!
def NM_separate_points_groups(hc, vc, instr, max_scanning_angle_range_deg = 4*10.1):

    ERRMSG_BEGIN = "Error separating NM measurement points: "
    global abs_meas_positions_ids
    n = len(hc)

    def grouping_not_succedded_result():
        # all_points_together = []
        # for i in range(n): all_points_together = all_points_together + [i]
        # group_indeces  = [all_points_together]
        return [ [ list(range(n)) ], ["all-points"], [ "o" ] ]

    NM_meas_part_symbols = [ 'o',  'x',   '.', 'v',  'o',  'x',  '.',  'v',  's' ]
    #nof_meas_parts_withouth_extra_points = 8   # 4x Declination series + 4x inclination series
    nof_scale_check_points = 0 if ((n % 8)  == 0) else 1
    #nof_extra_points = nof_scale_check_points
    group_indeces=[]
    extended_NM_meas_part_names = list(abs_meas_positions_ids)
    MAX_DIFF_FROM_FIXED_ANGLE = math.radians(0.5)  # Note: if acqusition time vc0 is different from the currently used vc0, then vc angles will not be very close to 90/270 degree.
    MAX_NR_OF_SCAN_POINTS = (n-nof_scale_check_points)  // 8
    MIN_NR_OF_SCAN_POINTS = 1  # max( (n-nof_scale_check_points), 1)

    assert MAX_NR_OF_SCAN_POINTS >= MIN_NR_OF_SCAN_POINTS
    act_meas_part_index = 0  # Index of decl / incl series
    meas_part_rotation_plane = "hc" # Supposing measurement started with a declination series
    # scan_start_index = 0   # Supporing that first point is a scan point of a declination series
    # meas_part_fixed_angle  = instr.get_horizontal_vc_angle("SU") if vc[0] < np.pi else instr.get_horizontal_vc_angle("SD")  # Telescope is in horizontal position (90 deg or 270 deg)
    nr_scan_points = -1 # Not yet known
    # assert ( vc[0] >= meas_part_fixed_angle - MAX_DIFF_FROM_FIXED_ANGLE ) and ( vc[0] <= meas_part_fixed_angle + MAX_DIFF_FROM_FIXED_ANGLE ), ERRMSG_BEGIN+"First measurement point is not on horizontal plane"
    extra_point_indeces = []
    min_nr_of_scan_point = MIN_NR_OF_SCAN_POINTS
    max_nr_of_scan_point = MAX_NR_OF_SCAN_POINTS

    i = 0
    while i < n-nof_scale_check_points:

        # We have to jump over the extra points and detect the start of the next scan cycle
        can_be_a_new_scan = (act_meas_part_index<8)  # Later we can change 'can_be_a_new_scan' only to "False"
        if (nr_scan_points>0) and (i > n-nof_scale_check_points-nr_scan_points):
            can_be_a_new_scan = False
        if meas_part_rotation_plane == "hc":
            # We are wtihin declination series
            if (abs(vc[i]-instr.get_horizontal_vc_angle("SU", "rad") ) > MAX_DIFF_FROM_FIXED_ANGLE)  and (abs(vc[i] - instr.get_horizontal_vc_angle("SD", "rad")) > MAX_DIFF_FROM_FIXED_ANGLE):
                # Declination series must be in horizontal plane
                can_be_a_new_scan = False
        if not can_be_a_new_scan:
            extra_point_indeces = extra_point_indeces + [ i ]
            i = i + 1
            continue
        # Testing if 'i' is a start of a scan
        scan_start_index = i
        # Checking later points
        if (scan_start_index == n - nof_scale_check_points - 1)  and ( nr_scan_points == 1)  :
            # Pure null method (no scanning)
            scan_end_index = n - nof_scale_check_points - 1
        else:
            for scan_end_index in range(scan_start_index+1, scan_start_index+1+max_nr_of_scan_point):
                hc_diff_deg = abs(angle_utils.angle_distance_deg( math.degrees( hc[scan_start_index] ), math.degrees( hc[scan_end_index] ) ) )
                vc_diff_deg = abs(angle_utils.angle_distance_deg( math.degrees( vc[scan_start_index] ), math.degrees( vc[scan_end_index] ) ) )
                if meas_part_rotation_plane == "hc":
                    if vc_diff_deg > math.degrees(MAX_DIFF_FROM_FIXED_ANGLE):
                        # Too big vc fixed angle change
                        scan_end_index = scan_end_index -1
                        break
                    if hc_diff_deg > max_scanning_angle_range_deg:
                        # Too big scan angle change
                        scan_end_index = scan_end_index -1
                        break
                else:
                    if hc_diff_deg > math.degrees(MAX_DIFF_FROM_FIXED_ANGLE):
                        # Too be hc fixed angle change
                        scan_end_index = scan_end_index -1
                        break
                    if vc_diff_deg > max_scanning_angle_range_deg:
                        # Too big scan angle change
                        scan_end_index = scan_end_index -1
                        break
                if scan_end_index >= n- nof_scale_check_points-1:
                    # Last possible scan point
                    break
                if scan_end_index >= scan_start_index + max_nr_of_scan_point - 1:
                    # Last possible scan point
                    break

        new_nr_of_scan_points = scan_end_index - scan_start_index + 1
        if new_nr_of_scan_points < min_nr_of_scan_point:
            # 'i' (that is 'scan_start_index') is not a start scan index
            extra_point_indeces = extra_point_indeces + [ i ]
            i = i + 1
            continue

        # 'i' (that is 'scan_start_index') is a start scan index
        if (new_nr_of_scan_points<5):
            pass
        assert new_nr_of_scan_points <= max_nr_of_scan_point, f"{new_nr_of_scan_points} {max_nr_of_scan_point}"
        if nr_scan_points < 0:
            # This was the detected first full-scan
            nr_scan_points = new_nr_of_scan_points
            # Scan point number is the same for all 8 measurement parts
            min_nr_of_scan_point = nr_scan_points
            max_nr_of_scan_point = nr_scan_points

        elif nr_scan_points != new_nr_of_scan_points:
            assert False, ERRMSG_BEGIN+"Scan point number changed from {} to {}".format(nr_scan_points, new_nr_of_scan_points)
        scan_index_group = list(  range(scan_start_index,scan_end_index+1) )
        group_indeces = group_indeces + [ scan_index_group ]
        act_meas_part_index += 1
        if act_meas_part_index==4: meas_part_rotation_plane = "vc" # Going from declination series to inclination series
        i = scan_end_index + 1

    if act_meas_part_index != 8:
        return grouping_not_succedded_result()

    assert len(group_indeces) == 8

    # Adding scale check point group
    if nof_scale_check_points > 0:
        assert nof_scale_check_points == 1 # The next line prepared only for 1 scale check point:
        group_indeces = group_indeces + [ [ n - nof_scale_check_points  ] ]
    else:
        # Removing scale check point:
        extended_NM_meas_part_names = extended_NM_meas_part_names[:-1]
        NM_meas_part_symbols = NM_meas_part_symbols[:-1]

    # Adding extra points group
    if len(extra_point_indeces) > 0:
        group_indeces += [ extra_point_indeces ]
        extended_NM_meas_part_names += [ "Extra-point{}".format("s" if len(extra_point_indeces) > 1 else "" ) ]
        NM_meas_part_symbols += [ '1' ]

    assert len(group_indeces) == len(extended_NM_meas_part_names)
    assert len(group_indeces) == len(NM_meas_part_symbols)

    return { 'group_indeces':group_indeces, 'meas_part_names': extended_NM_meas_part_names, 'meas_part_symbols': NM_meas_part_symbols }, nr_scan_points  # Format required by the 'meas_point_groups' parameter of the 'lmfit_with_randominitialvalues_execute()' function (actually, it will be the 'point_groups' parameter of 'plot_util.plot_xy()' function)


# def determine_DIM_measurement_positions(hc_NED, vc_theo, instr): # azimuth_def,
#     n = len(hc_NED)
#     assert  n == len(vc_theo)
    
#     count_in_positions = np.empty((8),dtype=int)
#     count_in_positions[:] = 0
#     position_indeces = np.empty((n),dtype=int)
#     for i in range( n ):
#         quarter = calc_NED_quarter_direction(hc_NED[i],vc_theo[i],instr)
#         sensor_pos = instr.get_sensor_UD_position(vc_theo[i])
#         pos_name = quarter + sensor_pos
#         pos_index = abs_meas_positions_ids.index(pos_name)
#         count_in_positions[pos_index] += 1
#         position_indeces[i] = pos_index

#     return position_indeces, count_in_positions



def determine_DIM_measurement_positions(mag, hc_NED,  vc_theo, instr, precalc_quarters = None, precalc_FLFR = None):
    N = len(hc_NED)
    assert  N == len(vc_theo)
    indeces_per_pos = {k: [] for k in abs_meas_positions_ids} # "dict.fromkeys()" is not working, as it shares one list object
    meas_nr_per_pos = dict.fromkeys(abs_meas_positions_ids, 0) 
    
    inclination_indeces = []

    for i in range( N ):
        quarter = calc_NED_cardinal_direction(hc_NED[i],vc_theo[i],instr)
        if precalc_quarters is not None: 
            if quarter != precalc_quarters[i]:
                assert quarter == precalc_quarters[i]
        face_position = instr.get_face_direction(vc_theo[i], "LR") 
        if precalc_FLFR is not None:
            assert face_position == precalc_FLFR[i]
        pos_id = quarter + face_position
        assert pos_id in abs_meas_positions_ids
        indeces_per_pos[pos_id] += [i]
        meas_nr_per_pos[pos_id] += 1
        if quarter in ["S", "N"]: 
            inclination_indeces += [ i ]

    # Trying to auto determine the scale check point_
    if (len(inclination_indeces)>4) and  ( (len(inclination_indeces) % 4) == 1 ):
        expected_meas_nr_per_pos = len(inclination_indeces) // 4  # For scale check position, we will expect one more sample
        for incl_pos_id in abs_meas_positions_ids[ abs_meas_positions_incl_slice ]:
            if meas_nr_per_pos[incl_pos_id] == expected_meas_nr_per_pos + 1:
                # This position may contain the scale check
                max_M_sample_index = np.argmax( np.abs(mag[ indeces_per_pos[incl_pos_id] ]) )
                scale_check_point_index = indeces_per_pos[incl_pos_id][max_M_sample_index]
                indeces_per_pos[SCALECHECK_POS_ID] += [ scale_check_point_index ]
                meas_nr_per_pos[SCALECHECK_POS_ID] += 1
                indeces_per_pos[incl_pos_id].remove(scale_check_point_index)
                meas_nr_per_pos[incl_pos_id] -= 1                
                break

    measurement_balanced = all(  meas_nr_per_pos[index] == meas_nr_per_pos[abs_meas_positions_ids[0]] for index in abs_meas_positions_ids[abs_meas_positions_declincl_slice] )

    return indeces_per_pos, measurement_balanced


# A rough way to remove levels in residuals, maybe mostly due to eccentricity errors
def NMposition_mean_removed_residuals_calc(raw_resids):
    if len(raw_resids)<16: return raw_resids
    
    N = len(raw_resids)
    nr_pf_parts = 8
    assert (N % nr_pf_parts)==0
    N_NMpos = N // nr_pf_parts
    res = np.empty(N)
    for ipart in range(nr_pf_parts):
        NMpos_indeces = slice(ipart*N_NMpos, ipart*N_NMpos+N_NMpos)
        NMpos_mean = np.mean( raw_resids[NMpos_indeces] ) 
        res[NMpos_indeces] =  raw_resids[NMpos_indeces] - NMpos_mean
    return res


#def fullmagfit_results_header_str_get():
#    return "samples random_step outliers | resid[nT] | I[°] D[°] | scale Moffset[nT] M_x[nT] M_y[nT] M_z[nT] | epsilon[arcmin] delta[arcsec] ortho_err[arcsec] " # "| hc_excenterr_ampl[']   hc_excenterr_phase[°]  vc_excenterr_ampl['] vc_excenterr_phase[°]  ae2_omega ae2_ampl['] ae2_phase[°]"

# -------------------------------------------------------------------
# --- Processing null method D series measurements in a simple way -- 
# -------------------------------------------------------------------
# This assumes that we measured in 0 nT positions.
# Sensor-Up is the "Face-left" reference direction.
# 'mag_meas_positions' is a structure array, containing declination part abolute measuremens data
# 'posindex' is a dictionary, defining measurement order
def NM_declination_part_basic_evaluation(mag_meas_positions, posindex, start_MIRA_UP,start_MIRA_DN, end_MIRA_UP, end_MIRA_DN, azimuth_mark_NED, H )    :
    
    # Averaging MIRA angles, with UP sensor pos as reference direction
    MIRA_DN_TO_UP_shift =  -np.pi if (start_MIRA_UP < start_MIRA_DN)  else +np.pi
    MIRA_mean_for_SU_pos = np.mean([start_MIRA_UP, start_MIRA_DN + MIRA_DN_TO_UP_shift, end_MIRA_UP, end_MIRA_DN + MIRA_DN_TO_UP_shift ]) 
    assert abs(MIRA_mean_for_SU_pos - start_MIRA_UP)<0.1 # Ensuring that averag directs to MIRA SU

    # Caclulating magnetic meridian in theo reference frame by averaging 0 nT positions. 
    # Average will be perpendicular to the measurement 0 nT direction, so will point to magnetic North or Southe.
    D_theo_mean = np.mean(mag_meas_positions["hc"] ) 
    # 'D_theo_mean' can direct either to magnetic North or to magnetic South 
    
    EU = mag_meas_positions['hc'][posindex["EU"]]
    ED = mag_meas_positions['hc'][posindex["ED"]]
    WU = mag_meas_positions['hc'][posindex["WU"]]
    WD = mag_meas_positions['hc'][posindex["WD"]]
    
    
    # If East-Up measured theo angle was higher than West-Up theo angle, then 'D_theo_mean' (that is a value just 
    # between the two) directs to North, so it is ok. 
    # In the opposite case, D_theo_mean must bo rotate by 180 deg:
    if EU < ED: D_theo_mean = (D_theo_mean + np.pi) % (2*np.pi)
        
    # Self-check with other solution:
    D_check_1, D_check_2 = amdp_base_formulas.Meridian( ED, WD, EU, WU  )
    
    if (abs(D_check_1 - D_theo_mean) > 1e-6) and (abs(D_check_2 - D_theo_mean) > 1e-6):
        assert False
    
    D_self_calc = azimuth_mark_NED  - MIRA_mean_for_SU_pos + D_theo_mean
    
    D_self_calc = (D_self_calc + 2*np.pi) % (np.pi * 2.0) 
    
 
    D_check = amdp_base_formulas.Declination(azimuth_mark_NED, np.mean((start_MIRA_UP, end_MIRA_UP)), np.mean((start_MIRA_DN, end_MIRA_DN)),  ED, WD, EU, WU ) 
    if abs(D_check-D_self_calc)  > 1e-6:
        assert False
        
    if D_self_calc > np.pi: D_self_calc -= (np.pi * 2.0)         
    
    # From CYG sheet:
    M_offset_declination_part = H * (EU - WU + ED - WD) / 4.0
    M_offset_check = amdp_base_formulas.Offset_D(ED, WD, EU, WU, H)
    if abs(M_offset_declination_part - M_offset_check ) > 1e-3:
        assert False
        
    delta = amdp_base_formulas.Delta(ED,WD,EU,WU)
    
    return D_self_calc, MIRA_mean_for_SU_pos, M_offset_declination_part, delta


# Callback of LM method.
def fullmagfit_diminstrument_based_DI_model(X, M_scale, M_offset, epsilon, delta,
                                            hc_ecc_phase, vc_ecc_phase,  hc_ecc_ampl,  vc_ecc_ampl,
                                            I, D, ortho_error,
                                            # ae2_ampl, ae2_phase_, ae2_omega,
                                            # ae3_ampl, ae3_phase_, ae3_omega, consider_ae3,
                                            model_extra_params):
    # neglect_ortho_error = model_extra_params[0]  # 'bool' is converted to 'double, so '1.0' means 'true'
    hc_theo_measured = X[0]
    vc_theo_measured = X[1]
    F = X[2]

    assert isinstance(model_extra_params,dict)
    varcorr_method = model_extra_params['varcorr_method']
    variometer_data   = model_extra_params['variometer_data']
    F_data = model_extra_params['F_data']
    M_DIM = model_extra_params['M_DIM']

    # F should be a time series:
    assert len(F) == len(hc_theo_measured)

    hc_theo = hc_theo_measured # + ae1_hc_ampl * np.cos ( hc_measured + phase__to_rad(ae1_hc_phase_,1.0) )  + ae2_ampl * np.cos ( hc_measured * ae2_omega +  phase__to_rad(ae2_phase_, ae2_omega) )
    vc_theo = vc_theo_measured # + ae1_vc_ampl * np.cos ( vc_measured + phase__to_rad(ae1_vc_phase_,1.0) )  + ae2_ampl * np.cos ( vc_measured * ae2_omega +  phase__to_rad(ae2_phase_, ae2_omega) )

    if not np.isnan(hc_ecc_ampl): hc_theo = hc_theo + hc_ecc_ampl * np.cos( hc_theo_measured + hc_ecc_phase )
    if not np.isnan(vc_ecc_ampl): vc_theo = vc_theo + vc_ecc_ampl * np.cos( vc_theo_measured + vc_ecc_phase )

    # if consider_ae3:
    #     hc =  hc + ae3_ampl * np.cos ( hc_measured * ae3_omega +  phase__to_rad(ae3_phase_, ae3_omega) )
    #     vc =  vc + ae3_ampl * np.cos ( vc_measured * ae3_omega +  phase__to_rad(ae3_phase_, ae3_omega) )

    hpa_n = dim_azimuth.default_azimuth_def.h_theo_to_ned(hc_theo)
    vpa_n = dim_azimuth.v_theo_to_ned(vc_theo)

    # print("F:{:.1f} hc={:.1f}° vc={:.1f}°  I={:.2f}° D={:.2f}°  epsilon:{:.1f}' delta:{:.1f}'  Moffs={:.1f} nT scale={:.1f} nT   ".format(F, math.degrees(hc_measured[0]), math.degrees(vc_measured[0]), math.degrees(I), math.degrees(D), rad_to_arcsec(epsilon), rad_to_arcsec(delta), Moffset, scale ))

    # determining_ID_relative_to_varcorr_data_ID_variation = len(X)>3
    if varcorr_method == "LM/ID-var": # determining_ID_relative_to_varcorr_data_ID_variation:
        I_var=X[3]
        D_var=X[4]
        res =  amdp_base_formulas.dim_instrument_model(F, vpa_n, hpa_n, I-I_var, D-D_var, epsilon, delta, M_offset, M_scale)
    elif varcorr_method == "LM/dyn":
        varcorr_mag_sensor = dim_instr.mag_sensor_model_class(M_offset = M_offset, M_scale = M_scale, epsilon = epsilon, delta = delta)
        M_corr = calc_absmeas_variation_correction_values(varcorr_method, hpa_n, vpa_n, varcorr_mag_sensor, F_data, variometer_data[gu.COMPNAME_INCL_ANGLE], variometer_data[gu.COMPNAME_DECL_ANGLE] )
        res = amdp_base_formulas.dim_instrument_model(F, vpa_n, hpa_n, I, D, epsilon, delta, M_offset, M_scale) - M_corr
    else:
        res = amdp_base_formulas.dim_instrument_model(F, vpa_n, hpa_n, I, D, epsilon, delta, M_offset, M_scale)
    if (M_DIM is not None):
        # Target is zero
        res = res  - M_DIM
    return  res




def calc_NED_cardinal_direction(hc_NED,vc_theo,instr):
   cardinal_index =  int ( (((hc_NED + math.pi/4) %  math.tau) / (math.pi / 2.0) ))
   if not instr.is_telescope_in_reference_direction(vc_theo):  
       cardinal_index = (cardinal_index + 2) % 4
   return TELESCOPE_CARDINAL_DIRECTIONS[cardinal_index]
   



# Finds the same sensor vertical position after the the rotates 180 degree in horizontal circle.
# Finds the position where sensor is on the same side (UP or DOWN) than the starting position
# 'vc' is in theo frame, and in radian.
def get_complement_vc_angle(vc):
    return (2*math.pi - vc + math.tau)  % math.tau


# Calculates the estimated nT values for absolute measurement positions (defined by hc_NED / vc_NED / epsilon / delta) from a variometer recording.
# 'F'/'I'/'D' are variometer+scalar-mag readings at the time of abolute measurement positions.
# 'F' can be a constant value.
# Pillar difference (if know) must be applied for passed 'F'/'I'/'D'data prior to calling this fucntion.
def calc_absmeas_variation_correction_values(varcorr_method, hc_NED, vc_NED, varcorr_mag_sensor, F, Incl, Decl):

    assert isinstance(Incl, (np.ndarray))
    assert isinstance(Decl, (np.ndarray))
    assert isinstance(F, ( float, np.ndarray))  # F can be a constant or an array
    assert len(hc_NED) == len(vc_NED)
    assert len(hc_NED) == len(Incl)
    assert len(hc_NED) == len(Decl)

    M_variometer_at_ADS_sensor_positions_at_meas_time = amdp_base_formulas.dim_instrument_model(F, vc_NED, hc_NED, Incl, Decl, varcorr_mag_sensor.epsilon(), varcorr_mag_sensor.delta(), varcorr_mag_sensor.offset, varcorr_mag_sensor.scale)

    if varcorr_method in [ "avg", "LM/dyn"]:
        M_variometer_at_ADS_sensor_positions_at_ref_time = amdp_base_formulas.dim_instrument_model(np.mean(F), vc_NED, hc_NED, np.mean(Incl), np.mean(Decl), varcorr_mag_sensor.epsilon(), varcorr_mag_sensor.delta(), varcorr_mag_sensor.offset, varcorr_mag_sensor.scale)
        # Removing the average mag value from the values, to get variation around 0 nT:
        M_variometer_correction = M_variometer_at_ADS_sensor_positions_at_ref_time - M_variometer_at_ADS_sensor_positions_at_meas_time
    elif varcorr_method == "T[0]":
        reference_time_index = 0
        M_variometer_at_ADS_sensor_positions_at_ref_time = amdp_base_formulas.dim_instrument_model(F[reference_time_index], vc_NED, hc_NED, Incl[reference_time_index], Decl[reference_time_index], varcorr_mag_sensor.epsilon(), varcorr_mag_sensor.delta(), varcorr_mag_sensor.offset, varcorr_mag_sensor.scale)
        M_variometer_correction = M_variometer_at_ADS_sensor_positions_at_ref_time - M_variometer_at_ADS_sensor_positions_at_meas_time
    elif varcorr_method in ["none", "LM/ID-var" ]:
        M_variometer_correction = np.zeros( len(hc_NED) )
    else:
        assert False, f"Unknown variation correction method: {varcorr_method}"

    return M_variometer_correction


# def other_plane(plane):
#     if plane == ROTATION_PLANE_VERTICAL: return ROTATION_PLANE_HORIZONTAL
#     if plane == ROTATION_PLANE_HORIZONTAL: return ROTATION_PLANE_VERTICAL
#     assert False


