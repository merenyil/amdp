"""
Description
@author: László Merényi, Domján Ádám
Author  = László Merényi (LM), Adam Domjan (AD)
Name    = trigonometric_fitting.py
Date    = 2021 April 30

Library for sinus/cosinus fittings and some sort of harmonic analysis.

"""


import numpy as np
import math
# # import numexpr as ne                  # Faster numpy calculation
# from scipy.optimize import curve_fit  # For trigonometric method
from lmfit import Model  # For trigonometric method
import matplotlib.pyplot as plt
import sys   # For sys.platform() command
import datetime
# import ads_utils as au
import misc_utils as mu

# ----------------------------------------------------------------------------------------------------
# Function originally from "absolute magnetic data processing - trigonometric (cosine/sine) method"
# lmerenyi: added omega
# ----------------------------------------------------------------------------------------------------
def cosine_model(angle, phase, amplitude, offset, omega):
    """
    Callback function used to fit one component cosine curve to measured/synthetic data

    The general model is:

    dim_sensor_x = amplitude * np.cos(angle - shift) + offset

    (scenario 1)
    vertical_angle_in_theo_frame := 90.0000 degree
    angle                        := horizontal_angle_in_theo_frame
    amplitude                     = (+1) * scale_x * |F| * cos(I)
    shift                         = D_theo - delta
    offset                        = offset_x

    (scenario 2)
    vertical_angle_in_theo_frame := 270.0000 degree
    angle                        := horizontal_angle_in_theo_frame
    amplitude                     = (-1) * scale_x * |F| * cos(I)
    shift                         = D_theo + delta
    offset                        = offset_x

    (scenario 3)
    horizontal_angle_in_theo_frame := D_theo
    angle                          := vertical_angle_in_theo_frame
    amplitude                       = (+1) * scale_x * |F| * cos(delta)
    shift                           = pi/2 + Incl_ned - epsilon
    offset                          = offset_x

    (scenario 4)
    horizontal_angle_in_theo_frame := D_theo + pi
    angle                          := vertical_angle_in_theo_frame
    amplitude                       = (-1) * scale_x * |F| * cos(delta)
    shift                           = pi/2 - Incl_ned - epsilon
    offset                          = offset_x

    param angle:     hori/vert angle of the theodolite resolved in theodolite frame
    param shift:     check the shift formula in between scenarios
    param amplitude: check the amplitude formula in between scenarios
    param offset:    offset of the X component of the magnetic sensor of DIM instrument
    return:          is the magnetic sensor X component of DIM
    """


    # return ne.evaluate('amplitude * np.cos(angle - phase) + offset')  # Quick calc
    # print("DEB: ampl={:.2f} offset={:.2f}   y={:2f}".format(amplitude, offset, amplitude * np.cos(angle[0] * omega - phase) + offset) )

    # Note " + phase" has to be used here, to have same phase that 'cosine_fit_linreg()' returns
    return amplitude * np.cos(angle * omega + phase) + offset

# Triangle wave model (similar to cosine model)
def triangle_model(X , phase, amplitude, offset, omega):

    angle = X
    n = len(angle)
    result = np.zeros(n)
    for i in range(0,n):
        angle_norm = mu.move_to_0_2pi_range(angle[i]*omega+phase)
        if ( angle_norm <= np.pi ):
            # Decreasing section of the triangle:
            result[i] = offset + amplitude * ( 1.0 - 2.0 * angle_norm/np.pi)
        else:
            # Increasing section of the triangle:
            result[i] = offset + amplitude * ( -1.0 + 2.0* (angle_norm-np.pi)/np.pi)
    return result


# addition of two different frequency sinus
def cosine_model_2(x, common_offset,
                   ampl_0, phase_0, omega_0,
                   ampl_1, phase_1, omega_1):
    return  common_offset + \
            ampl_0 * np.cos(x * omega_0 + phase_0) +  \
            ampl_1 * np.cos(x * omega_1 + phase_1)


# addition of three different frequency sinus
def cosine_model_3(x, common_offset,
                   ampl_0, phase_0, omega_0,
                   ampl_1, phase_1, omega_1,
                   ampl_2, phase_2, omega_2):
    return  common_offset + \
            ampl_0 * np.cos(x * omega_0 + phase_0) +  \
            ampl_1 * np.cos(x * omega_1 + phase_1) + \
            ampl_2 * np.cos(x * omega_2 + phase_2)



# ----------------------------------------------------------------------------------------------------
# Fits a measured data set to the 'fit_mode_str' type model, using nonlinear optimization.
# First returned variable is a string - in case of success it is empty string, in case of error
#                                       it is the error message
# If fit model is 'cos', and 'omega_guess' is np.nan, then omega is not optimized -> in this case it is better to use
#   'cosine_fit_linreg' function
# ToDo: optimize for redccted phase (between 0...1)
# ----------------------------------------------------------------------------------------------------
def gen_fit_nonlin_opt(   fit_mode_str,
                          y, x,
                          phase_guess, amplitude_guess, amplitude_min, amplitude_max,
                          offset_guess, offset_min = np.nan, offset_max  = np.nan,
                          omega_initial = 1.0, omega_min = np.nan, omega_max = np.nan,
                          print_results = True, draw_plot = True):

    if ( len(x) != len(y) ):
        return "len(x) != len(y)",np.nan,np.nan,np.nan,np.nan

    if (phase_guess < 0.0) or (phase_guess>=2*np.pi):
        return "phase-guess has wrong value", np.nan, np.nan, np.nan, np.nan

    # Initialize the model, with initial guesses
    best_phase      = np.nan
    best_ampl  = np.nan
    best_omega      = np.nan
    best_offset     = np.nan
    best_chisqr     = np.nan

    def func_fit_model_solution_callback(resid_std, esti_lmfit):
        nonlocal best_phase, best_ampl, best_offset, best_omega
        best_phase  = esti_lmfit[0]
        best_ampl   = esti_lmfit[1]
        best_offset = esti_lmfit[2]
        best_omega  = esti_lmfit[3]
        # Old version: "Offset={:.3f} ampl={:.3f} phase={:8.4f}° omega={:.2f} chi_sqr={:.2f}".format( best_offset, best_ampl,  math.degrees(best_phase),  best_omega, best_chisqr))
        return "resid-std:{:.2f}  offset:{:.2f} omega:{:.2f} ampl:{:.2f} phase:{:.3f}°   ".format(resid_std, best_offset, best_omega, best_ampl, math.degrees(best_phase) )


    # If initial phase guess is far from the result, then optimization may fail (happens very rarely)
    # We have to try at least with two phases!
    # And we run more random tries for each phase guesses - just to be sure to have the best solution

    for phase_guess_tries_index in [ 0, 1 ]:

        if fit_mode_str == 'cos':
            model = Model(cosine_model)
        elif fit_mode_str == 'triangle':
            model = Model(triangle_model)
        else:
            assert False, "Not supported fit model: {}".format(fit_mode_str)

        PHASE_MIN = 0
        PHASE_MAX = mu.pi_x_2 - 1e-8
        model.set_param_hint('phase',     value=phase_guess,     min = PHASE_MIN,             max = PHASE_MAX)
        model.set_param_hint('amplitude', value=amplitude_guess, min = amplitude_min, max = amplitude_max)
        if np.isnan(offset_min) or np.isnan(offset_max):
            model.set_param_hint('offset',    value=offset_guess)
        else:
            model.set_param_hint('offset',    value=offset_guess,    min = offset_min,    max = offset_max)

        if ( np.isnan(omega_min) ):
            model.set_param_hint('omega', value=omega_initial,   min = omega_initial - 0.001,         max = omega_initial + 0.001, vary = False)
        else:
            model.set_param_hint('omega', value=omega_initial,   min = omega_min,      max = omega_max)

        model_params = model.make_params()

        NOF_RANDOM_STEPS = 30
        if phase_guess_tries_index == 1:
            reference_params = best_results # Using first phase_guess as reference param
            reference_params_descr = "first phase try's result"
            b_plot_best_result_residuals = draw_plot
            first_run_best_stddev = best_resid_std
            first_run_best_chisqr =  best_chisqr
        else:
            reference_params = None
            reference_params_descr = ""
            b_plot_best_result_residuals = False # At the first phase_try, we do not yet plotting

        best_resid_std, best_results, params_at_the_bound_str, bestres_random_try_no, nof_removed_outliers, best_residuals, best_model_result  = \
             lmfit_with_randominitialvalues_execute(NOF_RANDOM_STEPS, model, model_params, y, x,
             reference_params = reference_params, reference_params_descr = reference_params_descr, reference_params_notused_indeces = [],
             get_solution_str_callback = func_fit_model_solution_callback, b_plot_best_result_residuals = b_plot_best_result_residuals,
             b_print_new_random_step_results = print_results, b_print_progress = True)

        # Program self check: in the second 'phase_guess_tries_index' run, the first 'phase_guess_tries_index' run's best results are also re-run, as we set 'reference_params = best_results'.
        # So in the second run it should not be possible to  get worst result, then in the first run:
        # De lehet, hogy a 'first_run best_model_result.best_chisqr'-t és 'best_model_result.chisqr'-t kéne összahaosníltani???
        if phase_guess_tries_index == 1: assert best_resid_std <= first_run_best_stddev

        # Model at best result:
        model_params_set_from_results( model_params, best_results )
        best_fitted_points =  model.eval(params=model_params, X=x)

        result_ok = (best_phase >= 0.0 ) and (best_phase < mu.pi_x_2)

        # checking result (sometime the result is not ok!):
        if not result_ok:
            assert False, "Atgondolni mit kell tenni"
            # Writing to log
            result_ok = False
            t = datetime.datetime.utcnow()
            t = t.strftime('%Y-%m-%d %H:%M:%S')
            f = open("warnings-{:}fit.log".format(fit_mode_str), "a")
            f.write("--------------- {} ------------------\n".format( t) )
            f.write("sin-fit result error: phase_guess_tries_index={} phase={} amplitude={} offset={}\n".format(phase_guess_tries_index, result_phase, result_ampl, result_offset ) )
            f.write(" phase_guess:{} amplitude_guess:{} amplitude_min:{} amplitude_max:{} offset_guess:{} offset_min:{} offset_max:{}".format(phase_guess, amplitude_guess, amplitude_min, amplitude_max,offset_guess, offset_min, offset_max))
            f.write(" angles: {}:\n".format(x))
            f.write(" mag: {}:\n".format(y))
            f.close()
            print("WARNING: fit result error!")

        # ToDo: test these: chisqr, redchi, aic, and bic.
        # SINFIT_BETTER_RESULT_QUALITY_MIN_RATIO = 1.0 + 1e-10
        # if result_ok and ( (nr_good_results == 0) or ( ( result_chisqr < best_chisqr) and  ( best_chisqr / result_chisqr  > SINFIT_BETTER_RESULT_QUALITY_MIN_RATIO )   )   ):
        #     if (nr_good_results >0):
        #         t = datetime.datetime.utcnow()
        #         t = t.strftime('%Y-%m-%d %H:%M:%S')
        #         f = open("warnings-sinfit.log", "a")
        #         f.write("--------------- {} ------------------\n".format( t ))
        #         f.write("Found better estimates: {} -> {}\n".format(best_chisqr, result_chisqr ))
        #         f.close()
        #     nr_good_results = nr_good_results + 1
        #     best_chisqr     = result_chisqr
        #     best_phase      = result_phase
        #     best_ampl  = result_ampl
        #     best_offset     = result_offset
        #     best_omega      = result_omega

        #     best_model_result = model_result
        #     best_residuals = result_residuals

        phase_guess = (phase_guess + 1.0 *np.pi ) % mu.pi_x_2
        if phase_guess <= PHASE_MIN: phase_guess =  PHASE_MIN  + 1e-4
        if phase_guess >= PHASE_MAX: phase_guess =  PHASE_MAX  - 1e-4

    # Very rarely, there is no result:
    # if nr_good_results == 0:
    #     f = open("warnings-{:}fit.log".format(fit_mode_str), "a")
    #     f.write("!!! No result with the LM method !!!\n")
    #     f.close()
    #     return "no result", np.nan, np.nan, np.nan, np.nan

    # returning phase, amplitude, offset

    if print_results:
        print(func_fit_model_solution_callback(best_resid_std, best_results))

    # Plotting best result:
    # Plotting only on Windows
    if draw_plot and sys.platform.startswith('win'):
        # Alternative way: plot 'best_fitted_points' and 'best_residuals'
        plt.plot()
        # best_model_result.plot() # Plots both fit results and residuals
        best_model_result.plot_fit()
        plt.show()  # or plt.draw()
        best_model_result.plot_residuals()
        plt.show()
        plt.pause(0.0001)


    return  "", best_phase, best_ampl, best_offset, best_residuals, best_chisqr, best_omega


# ----------------------------------------------------------------------------------------------------
# Fits a measured data set to the cosine modell, using linear regression.
# From https://www.mathworks.com/matlabcentral/fileexchange/41246-sine-function-fit.
# "y = offset + ampl * cos(phase+x)"
# This implementation does not optimize for omega (frequency) - it is assumed that the received angles are in
# real radian units.
#
# 'do_not_calc_xresid_near_peaks__percentage': calculation of x residuals close to sinus maximum and minimum is
#  not accurate, because noise in Y data (typically: magnetic data) introduce a big error in x residual (typically: angle measurement's residual).
#  Percantage means excluded percantage of the sinus amplitude (sinus peak).
#  If 'do_not_calc_xresid_near_peaks__percentage' is at least 100, then X resid is not calculated at all
def cosine_fit_linreg(y, x, result_for_sinus = False, omega_corr = np.nan, print_results = False , do_orig_plot = True,  do_resid_plot = False, do_plot_for_x_resid = False, do_not_calc_xresid_near_peaks__percentage = 10.0, fig = None, subplot_start_rowindex = 1):

    def normalize_angle_diff(angle_diff):
        if (angle_diff < -np.pi):
            return angle_diff  + 2*np.pi
        if (angle_diff > +np.pi):
            return 2*np.pi - angle_diff
        return angle_diff

    # x and y can be either list or np.ndarray
    # This conversion is not necessary:
    # if type(i_x) is np.ndarray:
    #     # Need to convert from ndarray to list
    #     x = i_x.tolist()
    # else:
    #     x = i_x
    #     assert type(x) is list

    # if type(i_y) is np.ndarray:
    #     # Need to convert from ndarray to list
    #     y = i_y.tolist()
    # else:
    #     y = i_y
    #     assert type(y) is list

    n = len(x)


    # Solving linear equation system in the form of 'y = offset + b1*cos(x) + b2*sin(x)'
    #  where 'offset'/'b1'/'b2' are the unknowns and 'cos(x)' and 'sin(x)' are the known 'x' values.
    # Matrix for the solver (two sub-formats, depending omega correction has to be used or not):
    if np.isnan(omega_corr):
        # No need to multiply 'x':
        x_measured = x
    else:
        x_measured = x * omega_corr
        # xm = np.zeros(n)
        # for i in range(0,n):
        #     #   x[i] = x[i] * omega_corr note: this can not be used, as the 'x' array will be changed at the caller also (parameter passed as reference)
        #    xm[i] = x[i] * omega_corr

    m = np.vstack ( [ np.ones(n), np.sin(x_measured), np.cos(x_measured) ] ).T

    solution,resid,rank,s = np.linalg.lstsq(m, y, rcond=None)  # Corresponds to '[ones(size(x)),sin(d*x),cos(d*x)]\m_meas' in Matlab
    # or:    solution = np.linalg.solve(m.T.dot(m), m.T.dot(y))
    offset = solution[0]

    # Result is: y = solution[0] + solution[1]*np.sin(x)  + solution[2]*np.cos(x)
    # print( y - (solution[0] + solution[1]*np.sin(x)  + solution[2]*np.cos(x)) )

    ampl   = np.sqrt(solution[1]*solution[1] + solution[2]*solution[2])
    # if solution[1]*solution[2] < 0.0:        ampl = -ampl


    # c = atan2(b1,b2);
    # 'atan2' is the 4 quadrant version of atan, which handles the entire interval, [0,2*pi]. "
    phase = math.atan2(solution[2],solution[1])

    # This is the result for the 'SIN' function!
    phase = mu.move_to_0_2pi_range(phase)

    # Calculating residuals in two different point of view:
    # residual for Y, and residual for the angle (X).



    # We have a fitted sinus function that describes ideal 'y' for each x_measured: 'y_fitted=f(x_measured)', where 'f()' is a sine function.
    # 'y_fitted = offset + ampl * sin(x_measured +phase)'

    y_fitted_array = offset + ampl * np.sin(x_measured + phase)
    y_resid_array = y_fitted_array - y


    # Now we calculate residual for 'x', that is more complicated.
    # The inverse of this sinus is 'arcsin': 'x_fitted=f(y_measured)', where 'f()' is an 'arcsin()' function plus an offset correction.
    # In 'arcsin()', for one 'y_measured' point, there can be more 'x_fitted' values!
    #
    # We use this function:
    # (y_measured - offset) / ampl = sin(x_fitted + phase)
    #  =>
    # x_fitted = arcsin( (y_measured - offset) / ampl ) - phase
    #
    # But, because of periodicity of 'sin()', we have to use this from:
    #
    # x_fitted = M * arcsin( mod(ampl_ratio+1,2) -1 ) - phase + C
    # Where
    #     ampl_ratio =  ((y_measured - offset) / ampl   (should be in the range of -1....+1)
    #     M =  +1, if "normalized" 'x_measured + phase' is within -pi/2....+pi/2
    #     M =  -1, if "normalized" 'x_measured + phase' is within  pi/2....3/2*pi
    #     C is estimated from x_measured:
    #     C = trunc( (x_measured + phase + pi/2) / pi ) * pi

    x_resid_array = np.zeros(n)
    arcsin_M = np.zeros(n)
    arcsin_C = np.zeros(n)


    for i in range(0,n):

        # Checking if we are close to sinus peaks: if yes, not calculating residual

        if abs(  math.sin(x_measured[i]+phase ) ) > 1 - do_not_calc_xresid_near_peaks__percentage / 100.0:
            x_resid_array[i]=np.nan
            continue

        ampl_ratio = (y[i]-offset)/ampl
        # # We maximize to 1.00/-1.00 avoid "math domain errror" on noisy data:
        # if ampl_ratio > 1.00:
        #     ampl_ratio = 1.00
        # elif ampl_ratio < -1.00:
        #     ampl_ratio = -1.00

        if abs(ampl_ratio)>1.0:           # measured Y has error, it is above/below of the sinus amplitute, we can not handle this case
            x_resid_array[i]=np.nan
        else:

            arcsin_step = ((x_measured[i] + phase) + np.pi/2) // np.pi  # 'X // pi' equals to "trunc(X/pi)'
            arcsin_M[i] = 1 if (arcsin_step % 2) == 0 else -1
            arcsin_C[i] = arcsin_step * np.pi

            x_fitted = arcsin_M[i] * np.arcsin( np.mod(ampl_ratio+1,2) -1 ) - phase + arcsin_C[i]

            x_resid_array[i] = x_fitted - x_measured[i]

            # Old version:
            # There are two solutions for 'acos()', we select the better one:
            # if abs (math.cos (x[i] + phase) ) > 0.999:
            #     # In this section, the angle residuls have very high deviation
            #     x_resid_array[i] = np.nan
            # else:
            #     acos_1 =  math.acos(ampl_ratio)
            #     acos_2 = -math.acos(ampl_ratio)  + 2 * np.pi
            #     angle_diff_1 = normalize_angle_diff( acos_1 - x[i] - phase )
            #     angle_diff_2 = normalize_angle_diff( acos_2 - x[i] - phase )
            #     if abs(angle_diff_1) < abs(angle_diff_2):
            #         x_resid_array[i] = angle_diff_1
            #     else:
            #         x_resid_array[i] = angle_diff_2


    if len(resid)==0:
        resid_std_1 = np.nan # Tapasztalat: ha nagyon nagy a meres-szam, akkor resid valamirt ures!
    else:
        resid_std_1 = resid[0]/n    # todo: atgondolni, mi ez! Neha jo, neha nem!

    # 2025.04-es tapasztalat szintetikus adatoknál: 'resid_std_1' nem jó!
    resid_std_2 = np.std(y_resid_array)

    # The above algorithm fits a sine wave, but we need cosine, that is different with phase only:
    if not result_for_sinus:
        phase =  phase - np.pi/2.0
        phase = mu.move_to_0_2pi_range(phase)


    if print_results:
        print("Offset={:.3f} ampl={:.3f} phase={:8.4f}° resid_std_1={:.2f} resid_std_2={:.2f}".format( offset, ampl,  math.degrees(phase),  resid_std_1, resid_std_2))

    if fig != None:  plot_nrows, plot_ncols = fig.axes[0].get_subplotspec().get_gridspec().get_geometry()

    resid_plot_symbol = [ 'o', '-' ] [ 1 ]

    if do_orig_plot:
        if fig != None:
            # Plotting on an existing figure, into a new row. The first row may contain more plots: the colus is on 'plot_ncols' variable!
            # We join this plot cells!
            assert (subplot_start_rowindex > 0) and (subplot_start_rowindex <= plot_nrows)
            plt.subplot(subplot_start_rowindex, plot_ncols, ((subplot_start_rowindex-1) * plot_ncols+1, subplot_start_rowindex * plot_ncols ) )
            subplot_start_rowindex += 1

        plt.plot( np.degrees(x), y, 'o', color='black');
        plt.plot( np.degrees(x), y_fitted_array, '-', color='red');
        plt.xlabel('X [deg]')
        plt.ylabel('Orig Y and fitted Y')
        # 2025.-04-08-ig ez volt: if fig != None: plt.show()
        if fig is None: plt.show()

    if do_resid_plot:
        if fig != None:
            # Plotting on an existing figure, into a new row. The first row may contain more plots: the colus is on 'plot_ncols' variable!
            assert (subplot_start_rowindex > 0) and (subplot_start_rowindex <= plot_nrows)
            plt.subplot(subplot_start_rowindex, plot_ncols, ((subplot_start_rowindex-1) * plot_ncols+1, subplot_start_rowindex * plot_ncols ) )
            subplot_start_rowindex += 1
        plt.plot( np.degrees(x), y_resid_array, resid_plot_symbol, color='black');
        plt.xlabel('X deg]')
        plt.ylabel('Y-resid')
        # 2025.-04-08-ig ez volt:  if fig != None: plt.show()
        if fig is None: plt.show()

    if do_plot_for_x_resid:
        if fig != None:
            # Plotting on an existing figure, into a new row. The first row may contain more plots: the colus is on 'plot_ncols' variable!
            assert (subplot_start_rowindex > 0) and (subplot_start_rowindex <= plot_nrows)
            plt.subplot(subplot_start_rowindex, plot_ncols, ((subplot_start_rowindex-1) * plot_ncols+1, subplot_start_rowindex * plot_ncols ) )
            subplot_start_rowindex += 1
        x_resid_array_arcsec = [math.degrees(i)*3600.0 for i in x_resid_array]
        plt.plot( np.degrees(x), x_resid_array_arcsec, resid_plot_symbol, color='gray');
        plt.xlabel('X [deg]')
        plt.ylabel('X-resid (arcsec)')
        if fig != None: plt.show()


    # Tapasztalatok: resid_std_1 inkább a teljes peak-to-peak (esetenként annál jóval nagyobb), 'resid_std_2' közelebb áll az stddev-hez!
    #

    # ToDo: return cos_function instead: result_cos_function = cos_function( offset = offset, ampl = ample, phase_ = rad_to_phase_(phase,1.0) )

    return phase, ampl, offset, y_resid_array, resid_std_1, resid_std_2, x_resid_array


# "opt_omega_list": List with True/False, depending on if omega has to be refined (optimized) or not
# 'omega_list' order is important:  it must start with the most significant signal. Unknown omegas must be set as np.nan.
def multicomp_cosine_fit(omega_list, opt_omega_list, Y, X, result_for_sinus, omega_min_list = [], omega_max_list = []):

    n = len(omega_list)

    assert n == len(opt_omega_list), "omega_list and opt_omega_list sizes are different!"
    assert (n==2) or (n == 3) , "This function recently supports only 2 and 3 component sinus fitting"

    common_offset_initial_guess = np.mean(Y)

    # These variables stores the initial guesses
    ampl_list   =  []
    phase_list  =  []

    if (len(omega_min_list)==0):
        # If omega range was not supplied, then trying to find optimized omega's limits (min/max)
        # Supposing that omega is increasing in 'omega_list' (this is not always true, but the problems this function was created for has always this nature: the
        # more dominant a signal compomnent, the lower its frequencys).

        for i in range(0,n):
            if not opt_omega_list[i]:
                omega_min = omega_list[i] - 0.001
                omega_max = omega_list[i] + 0.001
            else:
                if i < n - 1:
                    if opt_omega_list[i+1]: # Next omega is fixed, so probably this component is far from that
                        omega_max = omega_list[i+1]*0.8
                    else:
                        omega_max = omega_list[i+1]*0.999
                else:
                    omega_max = omega_list[i]*10.0

                if i == 0:
                    omega_min = omega_list[i]/10.0
                else:
                    if opt_omega_list[i-1]: # Next omega is fixed, so probably this component is far from that
                        omega_min = omega_list[i-1]*0.8
                    else:
                        omega_min = omega_list[i-1]*0.5

            omega_min_list = omega_min_list + [ omega_min ]
            omega_max_list = omega_max_list + [ omega_max ]

    assert( len(omega_min_list) == len(omega_max_list) )

    # Decomposing known frequemcies to have a starting value:
    remaining_Y = Y

    for i in range(0,n):
        print(" {:d}. wave component's stand-alone sin fit: ".format(i+1), end='' )

        if opt_omega_list[i]:
            # As we doesnt know omega, we can make an initial guess for sin parameters only using nonlinear opt:
            phase_guess = 0.0
            amplitude_guess = np.std(remaining_Y) * math.sqrt(2.0)

            amplitude_min = 0.0
            amplitude_max = amplitude_guess * 10
            offset_guess = 0.0
            print("starting-omega:{:.2f} omega-range:{:.2f}...{:.2f} ampl_max:{:.1f}".format( omega_list[i], omega_min_list[i], omega_max_list[i], amplitude_max ), end='' )

            err_str, phase, ampl, offset, remaining_Y, chisqr, omega = gen_fit_nonlin_opt ("cos", remaining_Y, X, phase_guess, amplitude_guess, amplitude_min, amplitude_max, offset_guess, np.nan, np.nan, omega_list[i], omega_min_list[i], omega_max_list[i] )
            assert err_str == "", "'cosine_fit_nonlin_opt()' returned error"
            # todo: remaining Y kiirasa,
            # 'cosine_fit_nonlin_opt()' és 'cosine_fit_linreg()' közös paraméterekkel kezelni, és mindegyiknél legyen egy 'verb_mode'


            phase_list  = phase_list  + [ phase ]
            ampl_list   = ampl_list   + [ ampl ]
            omega_list[i] = omega
            print(" from nonlin opt: ampl={:f} phase={:f} omega={:f}".format(ampl,phase,omega))

        else:
            # As we know omega, we can make an initial guess for sin parameters using analytical form:
            print( " fixed_omega={:.3f} ".format(omega_list[i]), end = '' )
            phase, ampl, offset, remaining_Y, resid_std_1, resid_std_2, resid_x = cosine_fit_linreg(remaining_Y, X, result_for_sinus,  omega_list[i], True, True )

            # kiplotolni mi van itt

            phase_list  = phase_list + [ phase ]
            ampl_list   = ampl_list  + [ ampl ]

        if (i==0): common_offset_initial_guess  = offset  # The first frequency is the most dominant frequency, offset is determined mainly by this freq.

    print("-------- Together opt. initial values ----------")
    print("Initial offset: ", common_offset_initial_guess)
    print("phase_list:" , phase_list)
    print("ampl_list:" , ampl_list)
    print("opt_omega_list:" , opt_omega_list)
    print("omega_list:"  , omega_list)
    print("omega_min_list:"  , omega_min_list)
    print("omega_max_list:"  , omega_max_list)


    # The 'magic': optimizing all componetns in one step:


    # If initial phase guess is far from the result, then optimization may fail (happens very rarely)
    for phase_guess_tries_index in [ 1,2 ]:
        if n == 3:
            model = Model(cosine_model_3)
        elif n == 2:
            model = Model(cosine_model_2)
        else:
            assert False, "Not supported n"

        model.set_param_hint("common_offset",  value=common_offset_initial_guess )  # ToDo: set min/max

        for i in range(0,n):

            model.set_param_hint("ampl_{:d}".format(i),    value=ampl_list[i] )  # ToDo: think on min/max. Maybe: min = ampl_list[i] / 10.0, max = ampl_list[i] * 10.0
            model.set_param_hint("phase_{:d}".format(i),   value=phase_list[i],   min = 0.0,             max = 2 * np.pi)

            if opt_omega_list[i]:
                model.set_param_hint("omega_{:d}".format(i),   value=omega_list[i],   min = omega_min_list[i],      max = omega_max_list[i])
            else:
                # Not optimizing omega, it is a fixed value:
                model.set_param_hint("omega_{:d}".format(i),   value=omega_list[i],   vary = False)


        model_params = model.make_params()

        result = model.fit(Y, params=model_params, x = X)
        # print( result.fit_report() )


        esti_lmfit =  np.array([*result.best_values.values()])

        common_offset = esti_lmfit[0]
        result_index = 1
        for i in range(0, n):
            ampl_list[i] =  esti_lmfit[result_index]
            result_index = result_index +1
            phase_list[i] =  esti_lmfit[result_index]
            result_index = result_index +1
            omega_list[i] = esti_lmfit[result_index]
            result_index = result_index +1

        result_chisqr = result.chisqr

        residuals = np.zeros( len(X) )
        for i in range(0, len(X) ):
            if n==3:
                residuals[i] = Y[i] - cosine_model_3( X[i], common_offset, ampl_list[0], phase_list[0], omega_list[0],
                                                                       ampl_list[1], phase_list[1], omega_list[1],
                                                                       ampl_list[2], phase_list[2], omega_list[2])
            elif n == 2:
                residuals[i] = Y[i] - cosine_model_2( X[i], common_offset, ampl_list[0], phase_list[0], omega_list[0],
                                                                       ampl_list[1], phase_list[1], omega_list[1])
            else:
                assert False

        resid_std = np.std(residuals)

        print("3-component sin fit results try {}:".format(phase_guess_tries_index))
        print("  common offset: ", common_offset)
        for i in range(0, n):
            print("  {:d}. comp: ampl={:f} phase={:6.2f}°  omega={:f}".format( i, ampl_list[i], math.degrees(phase_list[i]), omega_list[i] ) )
        print("  chisqr:{:f}  resid-stddev:{:f}".format(result_chisqr, resid_std) )

        # Trying with other phase:
        for i in range(0,len(phase_list)):
            phase_list[i] = phase_list[i] + np.pi
            if (phase_list[i]>=2.0*np.pi ):
                phase_list[i] = phase_list[i] - 2.0 * np.pi

    return common_offset,  phase_list, ampl_list, omega_list, result_chisqr, resid_std, residuals


# ----------------------------------------------------------------------------------------------------
# Interface for 'tf.cosine_fit_nonlin_opt' , where input data is in angle_meas 2D array.
# 'angle_component_index' is '0' or '1', depending on the rotation axis.
# ----------------------------------------------------------------------------------------------------
def cosine_fit_nonlin_opt_(mag_meas, angle_meas, angle_component_index,
                          phase_guess, amplitude_guess, amplitude_min, amplitude_max,
                          offset_guess, offset_min, offset_max):

    if ( len(mag_meas) != len(angle_meas) ):
        return "len(mag_meas) != len(angle_meas)",np.nan,np.nan,np.nan,np.nan

    # Extracting rotating angle values (we don't care about the other angle)
    rotating_angle_meas = np.empty([len(mag_meas), 1])

    for i in range( 0, len(mag_meas) ):
        rotating_angle_meas[i] = angle_meas[i, angle_component_index ]

    # HIBÁS / TODO: 'tf' sehol nincs importálva/definiálva ebben a fájlban, és 'cosine_fit_nonlin_opt'
    # nevű függvény sem létezik a projektben (a hasonló, valódi függvény neve 'gen_fit_nonlin_opt',
    # más paraméter-sorrenddel: annak első paramétere egy 'fit_mode_str', pl. "cos", ami itt hiányzik).
    # Emellett a fenti hibaág (5 érték) és az itteni sikeres visszatérés (7 érték, ha 'gen_fit_nonlin_opt'-ra
    # javítanánk) érték-száma sem egyezne. JAVÍTANDÓ, mielőtt bárki használni kezdi ezt a függvényt!
    return tf.cosine_fit_nonlin_opt(mag_meas, rotating_angle_meas,
                                 phase_guess, amplitude_guess, amplitude_min, amplitude_max,
                                 offset_guess, offset_min, offset_max,
                                 np.nan)


