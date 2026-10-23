# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 18:50:43 2023

@author: Laszlo Merenyi

dIdD magnetometer data processing

"""
import numpy as np
import math
import dat_file_reader
import matplotlib.pyplot as plt
import plot_utils
import dim_utils as du
import geomag_utils as gu
import os


# There are different naming conventions for the deflected components.
# The goal os this function is the find the currently used deflected component name
def get_deflected_component_name(searched_comp, component_list):
    assert searched_comp in ["D+","D-","I+","I-"]
    if searched_comp in component_list: return searched_comp
    searched_comp.replace('+','p')
    searched_comp.replace('-','m')
    if searched_comp in component_list: return searched_comp
    raise ValueError( f'DIDD component "{searched_comp}" not found within components "{component_list}"' )


def didd_process_deflected_values(timelabel, F, Dp, Dm, Ip, Im, i0, d0, calc_comps = "IDXYZ"): # ~didd_deflected_values_to_ID_XYZ~

  # ToDo: megcsinalni:
  # case DIDDCalcMode of
  #   F1: begin Fa:=thisF;Fb:=thisF;end;
  #   F1F2: begin Fa:=thisF;Fb:=nextF;end;
  #   MEAN: begin Fa:=( (thisF+nextF) /2.0 ); Fb:=Fa;end;
  #   INTERPOLATION: begin Fa:=thisF+(nextF-thisF)/5.0*2.0; Fb:=thisF+(nextF-thisF)/5.0*4.0;end;
  #   else exit;
  # end;

  # data^.conv.F:=Fa; // Fb volt a korábbi C-s valtozatban

  assert calc_comps in ["IDXYZ", "HDZ", "dI", "dD", "ID", "I", "D"]

  #def calc_error(errmsg, timelabel):
  #      raise ValueError( "DIDD calculation error: {} at time: {:}. Comps:{:.2f}/{:.2f}/{:.2f}/{:.2f}/{:.2f}".format(errmsg, mu.exceltime_2_datetimestr(timelabel), F,Dp,Dm,Ip,Im  ) )

  F_D = F
  F_I = F

  Ip2  = np.square(Ip)
  Im2  = np.square(Im)
  F_I2 = np.square(F_I)
  F_D2 = np.square(F_D)

  aI = np.sqrt((Ip2+Im2-2.0*F_I2)/2.0)

  denominator  = 4.0 * aI * F_I
  DENOMINATOR_MIN_VALUE = 1e-5   # Denominator can not be zero, or a close-to-zero value
  denominator [ np.abs(denominator) < DENOMINATOR_MIN_VALUE  ]  = np.nan
  
  # if np.abs(denominator[denominator_min_index]) < 1e-5: calc_error('aI*F=0', denominator_min_index)

  sin_dI = (Ip2-Im2)/denominator;
  #max_index = np.argmax(np.abs(sin_dI))
  #if sin_dI[max_index] > 1:
  #    calc_error('abs(sin_dI)>1', timelabel)
  sin_dI [ np.abs(sin_dI) > 1 ] = np.nan
  
  dI = np.arcsin(sin_dI)

  if calc_comps=="dI": return dI

  I = i0 + dI
  if calc_comps=="I": return I

  Dp2 = np.square(Dp)
  Dm2 = np.square(Dm)
  aD = np.sqrt((Dp2+Dm2-2.0*F_D2)/2.0)

  denominator = 4.0 * F_D * np.cos(I) * aD
  denominator [ np.abs(denominator) < DENOMINATOR_MIN_VALUE  ]  = np.nan
  # if np.min(np.abs(denominator)) < 1e-5: calc_error('F*cos(I)*aD=0', timelabel)
  
  sin_dD  = (Dp2-Dm2)/denominator
  sin_dD [ np.abs(sin_dD) > 1 ]  = np.nan

  #if np.max(np.abs(sin_dD)) > 1:
  #    calc_error('abs(sin_dD)({:}) >1  at i0={:.3f}'.format(np.max(np.abs(sin_dD)), i0), timelabel)
  dD = -np.arcsin(sin_dD)
  
  if calc_comps=="dD": return dD

  D  = d0 + dD
  if calc_comps=="D": return D

  if calc_comps=="ID": return I,D

  H = F_I*np.cos(I)
  Z = F_I*np.sin(I)
  
  if max(Z)>48000:
      assert False

  if calc_comps=="IDXYZ":
      X = H * np.cos(D)
      Y = H * np.sin(D)

      return I,D,X,Y,Z
  
  assert False
  

