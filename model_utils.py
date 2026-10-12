# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 20:45:18 2024

@author: László Merényi 

Recommended way to create a model:

    def model_function(X, opt_param_1, opt_param_2, ..., model_extra_params):
        <implementation of model function >

    import lmfit
    from lmfit import Model
    import model_utils

    [...]

    model = Model(model_function, independent_vars = [ 'X', 'model_extra_params' ] )
    [...]

    # Ver-A) Adding a parameter without unit information and without optional unit conversion:
    model.set_param_hint('obsdaqC_scale_const',   value = 1.0, min = SCALE_MIN, max = SCALE_MAX, vary = not zero_offset_only_test )
    model_params = model.make_params()

    # Ver-B) Adding a parameter with unit information and with optional unit conversion:
    model_params = lmfit.Parameters()
    model_params.add_many( lmfit.Parameter('M_scale',  value = M_scale,  min = mag_sensor.scale_min, max = mag_sensor.scale_max,  vary = self.optimize_M_scale, user_data= ("°", mu.RAD_2_DEG )) )

    [...]
    model_extra_params =  [ neglect_ortho_error ]
    modelfunc_X = [ ... ]

    best_resid_std, best_results, params_at_the_bound_names, bestres_random_try_no, bestres_nof_removed_outliers, best_res_residuals= model_utils.lmfit_with_randominitialvalues_execute(nr_of_random_steps, model, model_params, mag_measured, modelfunc_X, reference_params = None, reference_params_descr = "", reference_params_notused_indeces = [], get_solution_str_callback = None, b_plot_best_result_residuals = b_plot_best_result_residuals, b_print_new_random_step_results = b_print_new_random_step_results, model_extra_params =  model_extra_params, meas_point_groups=meas_point_groups)


"""

import math
import numpy as np
import random
import plot_utils
import copy
import matplotlib.pyplot as plt

import misc_utils as mu
import cli_utils as cu

conv_radian_to_deg    = ["°", 1.0/(2*np.pi)*360.0, "{:.4f}"]
conv_radian_to_arcsec = ["\"", 1.0/(2*np.pi)*360.0*3600.0, "{:.1f}"]


# returns "", "min" or "max', depending on the optimization results is at which of the two bounds
def check_param_if_within_range(value, range_min, range_max, factor = 1e-4):
    # if math.isclose(param_value, par.min) : # , rel_tol = (par.max - par.min)  / 1e6
    assert not np.isnan(range_min)
    assert not np.isnan(range_max)
    assert range_max > range_min
    param_range = range_max - range_min
    param_close_to_edge_limit = param_range * factor # 1e-10 sometime is too small. 2024.10: 1e-5 was too small also
    if value < range_min + param_close_to_edge_limit:
        # Result is close to range minimum:
        return "min"
    else:
        if (value <= range_min):
            assert (value > range_min)
    if value > range_max - param_close_to_edge_limit:
        return "max"
    else:
        if value >= range_max:
            assert (value < range_max)
    return ""


# This function check LM model results: checks if one or more parameter is at the lower or
# upper bound: if yes, then info about parameter is added to a resulting string (deprecated return mode).
# If any parameter is at the bound limits, then the function returns an empty lists.
# ~model_param_range_check_str~
def model_results_param_range_check(model_params, best_results):
    # result_str = ""
    params_at_range_names   = []
    params_at_range_indeces = []
    for param_index, (pname, par) in enumerate(model_params.items()):
        if par.vary:
            at_bound_str = check_param_if_within_range(best_results[param_index], par.min, par.max)
            if at_bound_str != "":
                params_at_range_indeces += [param_index]
                params_at_range_names   += [pname]

    return params_at_range_names, params_at_range_indeces



# If a fit returned 'results', then it can be reloaded to the model by this function
def model_params_set_from_results(model_params, results, exclude_indeces = [] ):
    for param_index,(pname, par) in enumerate(model_params.items()):
        if not (param_index in exclude_indeces):
            par.value  =  results[param_index]

def model_parameter_user_format_element(par, custom_name = None):
    param_name =  par.name if (custom_name is None) else  custom_name
    param_unit = ""
    format_mask = None
    mult_for_unit = 1.0
    if par.user_data is not None:
        if isinstance(par.user_data,(list,tuple)):
            # 'user_data' can store this information in a list/tuple: unit, conversion-to-unit, format-mask
            param_unit = par.user_data[0]
            if len(par.user_data)>1:
                mult_for_unit = par.user_data[1]
                if len(par.user_data)>2:
                    format_mask = par.user_data[2]
        else:
            # 'user_data' is just the unit
            assert isinstance(par.user_data,str)
            param_unit = par.user_data

    if format_mask is None:
        # Auto determination of format mask
        format_mask = mu.auto_format_mask( np.array( (par.min*mult_for_unit, par.max*mult_for_unit) ) , param_name )

    return param_name, param_unit, format_mask, mult_for_unit



# Formats a string for a new model result
def gen_get_solution_str_callback(parameters, resid_std, results):
    result = "{:.3f}".format(resid_std)
    par_index = -1
    for pname, par in parameters.items():
        par_index += 1
        if not par.vary: continue
        name,unit,format_mask,mult_for_unit = model_parameter_user_format_element(par)
        if unit != "": unit = " "+unit
        result = result + " "+pname+"="+format_mask.format(results[par_index] * mult_for_unit)+unit
    return result

def get_one_result_from_results_array(model_params, results, param_name):
    param_index = list(model_params.keys()).index(param_name)
    val = results[param_index]
    name,unit,format_mask,mult_for_unit = model_parameter_user_format_element(model_params[param_name])
    val = val * mult_for_unit
    return val

# General function to execute LM optimizations.
# 'reference_params' can be supplied to be a random-step starting paramter from the second step, and also the result (without optimization)
# is printed out (if 'reference_params_descr' is not None).
# 'reference_params' can be:
# - 'None': in that case we do not take into account it.
# - a best solution parameters from a previous run
# - a solution for that we are interesting the result of - e.g. "zero amplitude" case
#  'reference_params_descr' is a short string that describes the role of 'reference_params'
# 'reference_params_notused_indeces' are the indeces in the parameter list that are not constant.
# - 'auto_extend_range_param_indeces' : if these optimizing parameters have the best value close to one side of their limits, then
#      the range is shifted, so that it will be centered around the best value, then new LM tries are run.
# - 'meas_point_groups': if points can be separated into groups (e.g. different part of an absolute magnetic measurement), then residuals will be marked by the group information

# From 2025.03.19, 'model_params' is an input and output parameter. On output, it contains the optimized values!
def lmfit_with_randominitialvalues_execute(nr_of_random_tries, model, model_params, Y, iX, reference_params, reference_params_descr, reference_params_notused_indeces = [],  get_solution_str_callback = None, b_plot_best_result_residuals = False, b_print_new_random_step_results = True, b_print_progress = True, auto_extend_range_param_indeces = [], model_extra_params = None, meas_point_groups = None, timelabels = None, plot_vertical_separators = None, show_results_table = True):

    MIN_PERCENTAGE_FOR_SIGNIFICANT_CHANGE = 0.1
    nof_parameters = 0

    starting_parameters_copy = copy.deepcopy(model_params)

    # Checking if initial values are within bounds:
    for pname, par in model_params.items():
        nof_parameters += 1
        par_value_str = "\"{:3f}\"".format(par.value)
        if par.vary:
            # In case of 'vary' variables, starting value can not equal to a bound, as variable can stack to the bound
            if np.isnan( par.value ):
                mu.warning("\"{:}\" varying paramter has nan value!".format(pname))
                continue
            if par.min == -math.inf:
                mu.warning("\"{:}\" varying paramter has no min value!".format(pname))
            if par.max == +math.inf: mu.warning("\"{:}\" varying paramter has no max value!".format(pname))
            if (par.value == par.min) or (par.value == par.max):
                mu.warning( "Varying parameter \"{:}\" starting value {:} is at one of the bounds of range {:.3f}...{:.3f}. Value may stick to this value!".format(pname, par_value_str, par.min, par.max))
            assert (par.value >= par.min), "Model parameter \"{}\" starting value is less than bound min: {} < {}!".format(par.name, par_value_str, par.min)
            assert (par.value <= par.max), "Model parameter \"{}\" starting value is more than bound max: {} > {}!".format(par.name, par_value_str, par.max)
        else:
            # In case of constans variables, we can be just at the bound
            if not np.isnan(par.value): # Sometimes we use 'nan' for parameters that do not have effect
                assert (par.value >= par.min), "Model parameter \"{}\" starting value is less than bound min: {} < {}!".format(par.name, par_value_str, par.min)
                assert (par.value <= par.max), "Model parameter \"{}\" starting value is more than bound max: {} > {}!".format(par.name, par_value_str, par.max)
            assert( not ( (nof_parameters-1)  in auto_extend_range_param_indeces ) )

    #  Checking iX, if it has arrays (required for Python library)
    for x_param_comp in iX:
        # It seems, that independent variables must be numpy arrays
        assert isinstance(x_param_comp, np.ndarray)
        assert len(x_param_comp) == len(Y)


    best_results = np.empty(nof_parameters)
    rough_best_results = np.empty(nof_parameters)

    best_residuals = np.empty( len(Y)  )
    best_resid_std = np.nan
    resid_std = np.nan
    bestres_random_try_no = 0
    nof_removed_outliers = 0
    best_model_result = None

    first_run = True

    if b_print_progress: progress_print = cu.progress_print_class()

    if get_solution_str_callback is None: get_solution_str_callback = gen_get_solution_str_callback

    # ---------------- START OF NESTED FUNCTION -------------------
    def lmfit_with_actualvalues_execute(i_randomtry, progress_line_header, model_params):

        new_best_results = False
        nonlocal first_run
        nonlocal resid_std, best_resid_std, best_results, bestres_random_try_no, best_residuals, best_model_result
        # nonlocal params_at_the_bound_names
        if b_print_progress: progress_print.newline("{} prev:{:.4f} best:{:.4f}".format( progress_line_header, resid_std, best_resid_std   ) )

        if model_extra_params is None:
            result = model.fit(Y, params=model_params,  X=iX)
        else:
            result = model.fit(Y, params=model_params,  X=iX, model_extra_params = model_extra_params)

        results =  np.array([*result.best_values.values()])

        res_residuals = result.residual
        resid_std = np.std(res_residuals)

        # ToDo: use a statistical parameter from 'result' instead of 'resid_std', to determine best result

        assert( len(best_results) == len(results) )
        assert( nof_parameters == len(results) )

        if (first_run) or (resid_std < best_resid_std):
            new_best_results = True
            best_resid_std = resid_std
            best_residuals[:] = res_residuals[:]
            bestres_random_try_no = i_randomtry
            best_results[:] = results[:]
            #if b_return_model_result_copy:
                # ToDo: check if deepcopy is really needed to store and return the model results
            #    best_model_result = copy.deepcopy(result)
            if (first_run):
                rough_best_results[:] = results[:]
                having_significantly_different_result = True
            else:
                # Checking if we have a new solution that is significantly different from last 'rough_best_results',
                # that is worth to be printed
                having_significantly_different_result = False
                for j in range(nof_parameters):
                    if mu.change_percent_get(rough_best_results[j] , results[j]) > MIN_PERCENTAGE_FOR_SIGNIFICANT_CHANGE:
                        rough_best_results[:] = results[:]
                        having_significantly_different_result = True
                        break

            if having_significantly_different_result:
                new_best_solution_str = get_solution_str_callback(model_params, resid_std, results)
                assert type(new_best_solution_str) == str
                params_at_range_names, params_at_range_indeces = model_results_param_range_check(model_params, best_results)
                if len(params_at_range_names) > 0: # params_at_range_str !="":
                    new_best_solution_str += " | params at bounds: {}".format(params_at_range_names)

                # if printed_best_solution_str != new_best_solution_str: # Sometimes we have a slgithly better solution, with practically the same results - we do not print it again
                if b_print_new_random_step_results:
                    if b_print_progress: progress_print.insert_static_line("  rnd.try:{:03d}  {}".format(i_randomtry, new_best_solution_str))
                # printed_best_solution_str = new_best_solution_str

        # print("\n{:.3f} <> {:.3f}\n".format(resid_std, best_resid_std))
        first_run = False
        return new_best_results
    # ---------------- END OF NESTED FUNCTION -------------------


    for i_randomtry in range(0,nr_of_random_tries):

        new_best_results = lmfit_with_actualvalues_execute( i_randomtry,  "> rnd:{}/{}".format(i_randomtry+1,nr_of_random_tries),model_params)

        # Additionally step for cycling  parameters that are found to be best at or close to one side of the range:
        # Re-running the optimization by moving the range half-size right
        # if new_best_results and (len(auto_extend_range_param_indeces)>0):
        #     while True:
        #         new_result = False
        #         iparam = -1
        #         for pname, par in model_params.items():
        #             iparam += 1
        #             if iparam in auto_extend_range_param_indeces:
        #                 assert(par.vary)
        #                 par_range = (par.max - par.min)
        #                 close_to_range_tolerance = par_range / 1e6
        #                 if (best_results[iparam] < par.min + close_to_range_tolerance) or (best_results[iparam] > par.max - close_to_range_tolerance ) :
        #                     # Result is close to range minimum or maximum:
        #                     # Running a new LM cycle, where the range is shifted by half-window
        #                     orig_best_value = best_results[iparam]
        #                     model_params_1 = copy.deepcopy(model_params)
        #                     model_params_set_from_results(model_params_1, best_results, [] )
        #                     model_params_1[pname].min = par.min + par_range / 2
        #                     model_params_1[pname].max = par.max + par_range / 2

        #                     if (model_params_1[pname].value < model_params_1[pname].min):
        #                         # model_params_1[pname].value = model_params_1[pname].value + par_range
        #                         # Some more changes is required, as if we start the optimization in the local mininum, then it will not chagne
        #                         # model_params_1[pname].value = model_params_1[pname].value - par_range / 2000.0
        #                         model_params_1[pname].value = par.max + 1.0/2000.0# Old max is the center of the new, shifted window
        #                     else:
        #                         model_params_1[pname].value = model_params_1[pname].value + par_range / 2000.0
        #                     print('\nCYCLIC_PARAMETER_AT_BOUND: orig:{} new-try:{}\n'.format(best_results[iparam], model_params_1[pname].value))

        #                     lmfit_new_solution = lmfit_with_actualvalues_execute(i_randomtry, False, "> param-at-bound-test {}".format(par.name), model_params_1)
        #                     if lmfit_new_solution:
        #                         # New solution found!
        #                         # Shifting back the param range into the original range
        #                         if best_results[iparam] > par.max: best_results[iparam] = best_results[iparam] - par_range
        #                         assert (best_results[iparam] > par.min) and  (best_results[iparam] < par.max)
        #                         print("....");
        #                         input("New best result due to cyclic step: {} => {}  CLOSE:{}".format(orig_best_value, best_results[iparam], "Y" if math.isclose(best_results[iparam], orig_best_value) else "N" ));
        #                         print("....");
        #                         if not math.isclose(best_results[iparam], orig_best_value):
        #                             new_result = True
        #         # If we have more such parameters, and one of them has a new result, then the others have to be retested again
        #         if (len(auto_extend_range_param_indeces) == 1) or (not new_result): break

        # Trying with new starting values:
        if (i_randomtry == 0) and (reference_params is not None) :
            # After using the received default parameters, we use the caller supplied previous best results
            model_params_set_from_results(model_params,reference_params, reference_params_notused_indeces )
        else:
            # Applying random starting values
            for pname, par in model_params.items():
                if par.vary: par.value = par.min + random.random() * (par.max - par.min)

    # End of random tries.

    # Checking if best result is close to bounding limit, that can mean, that the limit is wrong, or there is no good solution
    params_at_the_bound_names, params_at_the_bound_indeces = model_results_param_range_check(model_params, best_results)

    # Checking if there is (are) paramter(s) at the edge, that range can be extended
    # Intersection of the two list of 'params_at_the_bound_indeces' and 'auto_extend_range_param_indeces'
    do_auto_range_extension_param_indeces = list( set(params_at_the_bound_indeces).intersection(auto_extend_range_param_indeces))

    if len(do_auto_range_extension_param_indeces)>0:
        # We re-run the optimization after shifting the range around the optimized value
        model_params_1 = copy.deepcopy(model_params)
        model_params_set_from_results(model_params_1, best_results, [] )
        iparam = -1
        for pname, par in model_params.items():
            iparam += 1
            if (iparam in do_auto_range_extension_param_indeces):
                assert(par.vary)
                par_range = par.max - par.min
                model_params_1[pname].min = par.value - par_range / 2
                model_params_1[pname].max = par.value + par_range / 2
                # print('\nAUTO_RANGE_EXTENSION: orig:{} new-try:{}\n'.format(best_results[iparam], model_params_1[pname].value))
        # Making again random tries:
        for i_randomtry in range(0,nr_of_random_tries):
            if lmfit_with_actualvalues_execute( i_randomtry,  "> rnd-for-auto-range-extension-params:{}/{}".format(i_randomtry+1,nr_of_random_tries,),model_params_1):
                print("NEW BETTER RESULT: {}".format(best_results))
            # Applying random starting values
            for pname, par in model_params_1.items():
                if par.vary: par.value = par.min + random.random() * (par.max - par.min)

    # We have the best random try's result in "best_results" array.
    if b_print_progress: progress_print.close()

    # For information, we make a simple run with the user supplied 'reference_params', without optimization
    if (reference_params is not None) and (reference_params_descr is not None):
        param_index=-1
        for pname, par in model_params.items():
            param_index += 1
            if not (param_index in reference_params_notused_indeces):
                par.value  =  reference_params[param_index]
            else:
                par.value  =  best_results[param_index]
        modelfunc_res_at_reference =  model.eval(params=model_params, X=iX)
        assert ( len(modelfunc_res_at_reference) == len(Y)  )
        residual_of_reference = np.empty(len(Y))
        residual_of_reference = modelfunc_res_at_reference - Y
        # for i in range(0,len(Y)): residual[i] =  modelfunc_res_at_prevres[i] - Y[i]
        resid_std = np.std(residual_of_reference)
        if b_plot_best_result_residuals:
            print("  {}-res:  {}".format(reference_params_descr, get_solution_str_callback(model_params, resid_std, reference_params) ))
            yaxis_label = 'residual at ref. parameters (\"{}\")'.format(reference_params_descr)
            plot_utils.plot_xy(range(len(residual_of_reference)), residual_of_reference, 'LMFIT - sample number', yaxis_label, 'x' );


    if b_plot_best_result_residuals:
        nr_of_resid_plots = 1 if timelabels is None else 2

        plt.figure( figsize=(10, 4) )
        plt.suptitle('Residual of the best solution')
        plt.subplot(nr_of_resid_plots,1,1)
        # Plotting residuals in simeple increasing order
        plot_utils.plot_xy(range(len(best_residuals)), best_residuals, 'LMFIT sample number', 'residual', 'x', point_groups = meas_point_groups, close_plot = False, vertical_separators =plot_vertical_separators )
        # Plotting residuals over time (if timelabels are supplied):
        if timelabels is not None:
            assert len(timelabels) == len(best_residuals)
            plt.subplot(2,1,2)
            plot_utils.plot_datetime_y( timelabels, [ best_residuals], "residuals" , 'x',  b_close_plot = False  )
        plt.show()


    # Creating a table with the result and starting values.
    # A much more simple but less detailed method is use: 'model_params.pretty_print()'
    if show_results_table:
        ttable = cu.tabulate_table_creator()
        for param_index,(pname, starting_par) in enumerate(starting_parameters_copy.items()):
            ttable.add_from_model_params_(starting_par,  optimized_value = best_results[param_index], starting_value = starting_par.value )
            model_params[pname].value = best_results[param_index]  # This does the same as 'model_params_set_from_results()'
        ttable.print()

    return best_resid_std, best_results, params_at_the_bound_names, bestres_random_try_no, nof_removed_outliers, best_residuals

