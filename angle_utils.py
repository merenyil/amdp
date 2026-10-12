# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 17:22:49 2025

@author: Merényi László
"""

import numpy as np
import math

angle_expression_decimal_digits = {"deg": 4}
user_angle_unit = "deg"
user_angle_unit_decimal_digits = angle_expression_decimal_digits[user_angle_unit]

angle_comparison_tolerance = { "rad": 1e-5 }

RAD_2_ARCMIN = 360.0*60.0 / math.tau
RAD_2_ARCSEC = 360.0*60.0*60.0 / math.tau


def rad2gon(rad):
    return rad / np.pi * 200.0

def gon2rad(gon):
    return gon / 200.0 * np.pi

def rad2arcsec(rad):
    return rad / np.pi * 180.0 * 3600.0

def rad_to_arcsec(rad):
    return np.degrees(rad)*3600.0

def arcsec_to_rad(arcsec):
    return np.radians(arcsec/3600.0)

def rad_to_arcmin(rad):
    return np.degrees(rad)*60.0

def arcmin_to_rad(min):
    return np.radians( min  / 60.0 )

def angle_distance_deg(a1_deg, a2_deg):
    a_deg = a2_deg - a1_deg
    a_deg = (a_deg + 180) % 360 - 180
    return a_deg

def angle_distance(a1, a2):
    a = a2 - a1
    a = (a + math.pi) % math.tau - math.pi
    return a

def select_equivalent_angle(angle, reference_angle, max_error = None, assert_if_error_limit_exceeds = True):
    #if isinstance(angle, np.ndarray):
    #    result = np.copy(angle)
    #    for i in range(result.shape[0]):
    #        result[i] = select_equivalent_angle(result[i], reference_angle, max_error, assert_if_error_limit_exceeds)
    #    return result

    # while angle > reference_angle + np.pi: angle -= np.pi*2
    # while angle < reference_angle - np.pi: angle += np.pi*2
    
    angle = (angle - reference_angle + math.pi) % math.tau + reference_angle - math.pi
    
    if max_error is not None:
        angle_err = np.max(abs(angle_distance(angle,reference_angle)))
        if angle_err>max_error:
            is_array = isinstance(angle, (np.ndarray, list, tuple) )
            if assert_if_error_limit_exceeds:
                errmsg = "Equivalent angle calculation error"
                if not is_array: 
                    errmsg += ": {:.1f} / {:.1f} err:{:.3f} limit:{:.3f}".format(np.degrees(angle), np.degrees(reference_angle), np.degrees(angle_err), np.degrees(max_error) )
                raise ValueError(errmsg)
            if not is_array:
                angle = np.nan
            else:
                angle = np.full_like(angle, np.nan, dtype=float)
            
    return angle


def angle_avg(a1,a2):
    # Called for DIM ref. mark readings angles. 'a1' and 'a2' are transferred to the same direction, so they are close to each other.
    # Eventual 360 deg wrap around is handled!
    # Simple equation without wrap-around handling: (a1+a2)/2.0
    # https://stackoverflow.com/questions/491738/how-do-you-calculate-the-average-of-a-set-of-circular-data
    res =  (np.arctan2( np.sin(a1) + np.sin(a2), np.cos(a1) + np.cos(a2) ) + math.tau) % math.tau
    if np.any(abs(a1-a2)<0.1):
        # Typical, not wrap around situation
        diff_check = abs(res - (a1+a2)/2.0)
        assert np.any(np.degrees(diff_check) < 0.001)
    return res


# Converts a string to radian.
# If the string can be a single number, then it is assumed to be originally in radians
# The string can have these units:
#  "deg"
#  "rad"    
# ToDo: extend the list
def str_to_radians(value_str):
    
  value_str = value_str.strip()
  if value_str.lower().endswith("deg"):
       number_part = value_str[:-3].strip()  
       degrees = float(number_part)
       return math.radians(degrees)
  if value_str.lower().endswith("rad"):   
      value_str = value_str[:-3].strip()  
  return float(value_str)    
    

def user_angle(rad):
    assert user_angle_unit == "deg"
    return np.degrees(rad)


def add_two_cosine_with_same_omega(ampl1, phase1, ampl2, phase2):
	# From: https://dspguru.com/files/Sum_of_Two_Sinusoids.pdf
	ampl_sum = np.sqrt( pow(ampl1*np.cos(phase1) + ampl2*np.cos(phase2),2)  + pow(ampl1*np.sin(phase1) + ampl2*np.sin(phase2),2) )
	phase_sum = np.arctan(  (ampl1*np.sin(phase1)+ampl2*np.sin(phase2))  / (ampl1*np.cos(phase1) + ampl2*np.cos(phase2)) )
	return ampl_sum, phase_sum
