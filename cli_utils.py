# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 21:13:22 2026

@author: Laszlo Merenyi

Function for decorated outputting information for the user, and
so user input
"""

import tabulate
import model_utils # Required for 'tabulate_table_creator'
import sys

wait_for_enter_at_warnings = True    # Gloval for the project, can be modified from outside
wait_for_enter_at_infos = True       # Gloval for the project, can be modified from outside 


def wait_for_enter(force = False):
    if wait_for_enter_at_infos or force:
        input("Press Enter to continue...")

def warning(msg, allow_skipping_repeated_warnings = False, serious = False):

    if not isinstance(allow_skipping_repeated_warnings, bool):
        assert isinstance(allow_skipping_repeated_warnings, int)
        # 'allow_skipping_repeated_warnings' is an integer variable: first 'allow_skipping_repeated_warnings' characters of warning message is constant.
        msg_compare_part = msg[:allow_skipping_repeated_warnings]
        allow_skipping_repeated_warnings = True
    else:
        # 'allow_skipping_repeated_warnings' is booolean
        assert  isinstance(allow_skipping_repeated_warnings, bool)
        msg_compare_part =  msg

    print("Warning: ",msg)
    if allow_skipping_repeated_warnings and (msg_compare_part in warning.skipped_warning_messges): return
    if (not wait_for_enter_at_warnings) and (not serious): return

    if allow_skipping_repeated_warnings:
        while True:
            resp = input("Press [Enter] to continue, or type [s]+[Enter] to skip this warnings")
            if resp == "" : return
            if resp in [ "s", "S" ]:
                 warning.skipped_warning_messges += [msg_compare_part]
                 return
    else:
        input("Press [Enter] to continue...")

warning.skipped_warning_messges = []


def error(msg, is_fatal = True):
    print("ERROR: {}".format(msg))
    if is_fatal:
        sys.exit()
    else:
        input("Press [Enter] to continue...")
        
        
class progress_print_class:

    def __init__(self):
        self.last_line_len = 0

    def init(self):
        self.last_line_len = 0


    def delete_last_line(self):
        if self.last_line_len>0:
            # print("\r", end='')
		    # delete_spaces_after = self.last_line_len - len(line)
            print('\b' * self.last_line_len, end = '', flush=True)
            print(' '  * self.last_line_len, end = '', flush=True)
            print('\b' * self.last_line_len, end = '', flush=True)
            print('\r', end='')
            self.last_line_len = 0

    def newline(self,line):
        self.delete_last_line()
        # else:
        #    delete_spaces_after = 0
        line = '['+line+']'
        print (line, end='', flush=True)
        # delete_spaces_after>0: print(" "*delete_spaces_after, end='', flush = True)
        self.last_line_len = len(line)

    def close(self):
        self.delete_last_line()
    	# if (self.last_line_len>=0):
    	# 	print ("\r", end='')
    	# 	print(" "*self.last_line_len, end='')
    	# 	self.last_line_len = -1
    	# 	print("\r", end="")

    def insert_static_line(self, line):
    	# if self.last_line_len>=0:
    	# 	print ("\r", end='')
    	# 	delete_spaces_after = self.last_line_len - len(line)
    	# else:
    	# 	delete_spaces_after = 0
        self.delete_last_line()
    	# print(line, end = '')
    	# if  delete_spaces_after>0: print(" "*delete_spaces_after, end='', flush = True)
    	# print("")
    	# self.last_line_len = -1
        print(line, flush=True)


# 'possible_responses is a list, where each element is list, with two elements:
# response and message
def ask_user(msg, possible_responses, allow_simple_enter = False):

    accepted_responses = []
    for resp in possible_responses:
        assert len(resp)==2
        msg=msg+" {}:{}".format(resp[0], resp[1])
        accepted_responses = accepted_responses + [resp[0]]
    msg = msg + " >"


    while True:
        resp = input(msg)
        if resp in accepted_responses:
            return resp
        if allow_simple_enter and resp=="":
            return resp


# This utility class is used to add a record in a list-of-list, that is used for print a 'tabulate' table
# It can do the print directly.
class tabulate_table_creator:
    def __init__(self, headers = ["Param.name", "unit", "Value", "initial\nval", "range\nmin", "range\nmax", "out of\nrange"]):
        self.params_at_the_bound_names = []
        self.result_tabulate_data = []
        self.headers = headers

    def add(self, param_name, param_unit, format_mask, value, mult_for_unit = 1, starting_val = None, range_min_val = None, range_max_val = None, check_range_out = True):
        range_out_str = ""
        if starting_val is None:
            starting_val_str = ""
        else:
            starting_val_str = format_mask.format(starting_val * mult_for_unit)
        assert (range_min_val is None)  == (range_max_val is None)
        if range_min_val is None:
            range_min_str = ""
            range_max_str = ""
        else:
            range_min_str = format_mask.format(range_min_val*mult_for_unit)
            range_max_str = format_mask.format(range_max_val*mult_for_unit)
            if check_range_out:
                range_out_str = model_utils.check_param_if_within_range( value*mult_for_unit, range_min_val*mult_for_unit, range_max_val * mult_for_unit )
                if range_out_str != "":
                    self.params_at_the_bound_names += [ param_name ]
        self.result_tabulate_data += [ [ param_name, param_unit, format_mask.format(value*mult_for_unit), starting_val_str, range_min_str, range_max_str, range_out_str ] ]

    # If 'value' is nan, then value is read from 'model_param_item.value', and it that case,
    # 'mu.model_params_set_from_results(model_params, best_results)' must be called!
    # If 'format_mask' is None, then it is automatically determined
    def add_from_model_params_(self, par, optimized_value, starting_value = None, param_unit = "", format_mask = None, rename = None ):
        param_value = optimized_value if optimized_value is not None else par.value
        param_name,param_unit,format_mask,mult_for_unit = model_utils.model_parameter_user_format_element(par, rename)
        if par.vary: param_name +=" *"
        self.add(param_name, param_unit, format_mask, param_value, starting_val = starting_value, range_min_val = par.min, range_max_val = par.max, mult_for_unit = mult_for_unit)

    def print(self, table_type = 'pretty'):
        print(tabulate.tabulate(self.result_tabulate_data, headers= self.headers,  tablefmt=table_type))  # 'simple' or 'fancy_grid' or 'pretty'

# Some simiar functions:
    
def print_text_in_frame(main_line,second_line=None,border = "---"):
    main_line =border+" "+main_line.strip()+" "+border
    print(border[0] * len(main_line))
    print(main_line)
    if second_line is not None:
        second_line = border+" "+second_line+" "
        second_line = second_line + border[0]*(len(main_line)-len(second_line))
        print(second_line)
    print(border[0] * len(main_line))


def print_framed(text, ch='*', pad=1):
    lines = text.split('\n')
    width = max(len(l) for l in lines) + 2 * pad
    border = ch * (width + 2)
    print(border)
    for l in lines:
        print(f"{ch}{' ' * pad}{l.ljust(width - 2 * pad)}{' ' * pad}{ch}")
    print(border)     
    
    
def print_header(s):
    linelen=max(80,len(s)+4)
    before_text = (linelen - len(s) - 2) // 2
    after_text   = linelen - len(s) - 2 - before_text
    
    print('-'*linelen)
    print('-'*before_text,' ',s,' ','-'*after_text)
    print('-'*linelen)
    