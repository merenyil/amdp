# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 20:56:41 2023

@author: Merényi László, Mingeo Ltd.

Python library to read in maglin saved compressed binary files.
Main function is 'sc_load()'.
Please see notes at this function.


Copyright (c) 2023-2024 Mingeo Ltd., László Merényi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import os
import datetime
import numpy as np


MAX_STRING_HEADERFIELD_LENGTH   = 255
MAX_FILE_FORMAT_VERSION_LEN     = 10
MAX_HEADER_PARAMETER_NAMES_LEN  = 120  # For station/instrument and DAQ names
RECORD_BUFFER_INCREASE_STEP     = 65536  # To make faster memory allocation, we increase data buffer size not one-by-one, but with higher steps
CHUNK_SIZE = 65536    # Processed input data file is read by this number of bytes
EXCEL_DATETIME_EPOCH = datetime.datetime(1899, 12, 30)    # Note: not "31st Dec but 30th"!


class component_class:
    name = ""
    unit = ""
    compr_base = None # it becomes an integer number after the first data record
    number_total_width = 0
    nr_of_decimals = 0
    mult_factor_for_decimals = 1
    stat_count = [0,0,0,0,0]


# 'scload' is the main function of this unit: it reads in a compressed binary file into memory.
#
# Input parameters
# ----------------
# - 'file_name': full file name
# - 'add_msec_to_data': a boolean flag defining is millisecond part of the time label shoud be
#    reported as a component
# - 'comp_indeces': indexes of the required components. If it is None, then all components from the
#   file are read and returned. If it is a list or tuple, then only the component with the listed
#   indeces (started from '0') are stored and returned.
# In case of sub-second file: the first component (index 0) is the component after the "msec" pseudo-component.
#
# Returned values
# ---------------
#  'errmsg, warnings, timelabels, records, stored_comps'
# - 'errmsg': 'None' if there was not error, or a string, if there was an error
# - 'warning': list of warning messages, or empty list
# - 'timelabels': numerical representation of time labels, in a a "n x 1" shaped array, where
#     'n' is the record number. Spreadsheet programs generally recognized these numbers as date-time values.
# - 'records': data in an "n x c" shaped array, where 'c' is the number of components.
# - 'stored_comps': a list of 'component_class' instances, defining the loaded components.
def sc_load(file_name, add_msec_to_data = True, comp_indeces = None):

    ERRMSG_UNEXPECTED_EOF = "Unexpected EOF"

    f = None
    chunk=[]
    chunk_len = 0
    index_within_chunk=0

    nr_of_processed_sections = 0
    verblevel = 0
    warnings = []
    nr_of_total_data_records = 0
    nr_of_section_data_records = None
    nr_of__0x20_second_time_jump = 0

    # Variables for time label calculation
    msec_in_first_datacolumn=None
    do_second_step_due_to_msec_overflow = False
    record_time = None  # Stores timelabels in one-second resolution
    record_time_msec = None  # Comes from an additional component, added to record_time in case of subsecond files
    default_next_msec = 0
    default_time_step__sec = None
    default_time_step__msec = None


    act_record_index = -1
    timelabels = None
    records = None

    # Parameters coming from section headers:
    file_format_version = None
    station_name = None
    instrument_name = None
    DAQ_app_name = None
    # time_step__timedelta = None
    all_comps = []
    nr_of_all_comps = None
    stored_comp_indeces = None # A modified version of user supplied 'comp_indeces'. Indeces in this list corresponds to files original component indeces ('0' means millisecond in sub-second files)

    try:

        def get_next_byte():
            nonlocal f, index_within_chunk, chunk, chunk_len
            global CHUNK_SIZE
            if index_within_chunk<chunk_len:
                index_within_chunk += 1
                return int(chunk[index_within_chunk-1])
            chunk = f.read(CHUNK_SIZE)
            if not chunk:
                raise Exception( ERRMSG_UNEXPECTED_EOF )
            index_within_chunk = 1
            chunk_len = len(chunk)
            return int(chunk[0])

        def read_section_params():
            nonlocal nr_of_processed_sections, file_format_version, warnings, station_name, instrument_name, DAQ_app_name, default_time_step__sec, default_time_step__msec, msec_in_first_datacolumn, record_time, default_next_msec
            nonlocal nr_of_all_comps, all_comps
            nonlocal do_second_step_due_to_msec_overflow
            nonlocal timelabels, records, act_record_index
            nonlocal add_msec_to_data, comp_indeces, stored_comp_indeces

            section_header_errmsg = "section header error: "

            # Reads a string type section parameter
            def read_string(max_char_count = MAX_STRING_HEADERFIELD_LENGTH):
                result = ""
                for i in range(max_char_count + 1):
                    b = get_next_byte()
                    if b == 0:
                        if result == "": raise Exception( section_header_errmsg+"empty section variable")
                        return result
                    result = result + chr(b)
                # The string should already be terminated
                raise Exception( section_header_errmsg+"too long string variable in section header")

            expect_compr_base_values = (nr_of_processed_sections==0)
            nr_of_processed_sections += 1
            if (verblevel>0): print(" ** Section no. %{} **".format(nr_of_processed_sections))
            # Section should begin by 4 x uint8_t(0)
            for i in range(4):
                if get_next_byte() != 0:
                    raise Exception(section_header_errmsg+"#1")

            while (1):

                section_param_chr = chr(get_next_byte())

                if section_param_chr == 'v':
                    new_file_format_version = read_string(MAX_FILE_FORMAT_VERSION_LEN)
                    if file_format_version is None:
                        if verblevel>0: print("  File format version: ""{}""".format(new_file_format_version))
                    elif file_format_version != new_file_format_version:
                        warnings = warnings + ["File format version changed within file: ""{}"" ->""{}""".format(file_format_version, new_file_format_version) ]
                    file_format_version = new_file_format_version
                    if file_format_version> "SC2.1":
                        raise Exception("File format version {} is not supported".format(file_format_version)  )
                    continue

                if section_param_chr == '>':  # ">>": data section begins...
                    if chr(get_next_byte()) != ">" : raise Exception(section_header_errmsg+"#22")
                    if record_time is None: raise Exception(section_header_errmsg+"missing date field")
                    break # End of section parameters

                # Parsing section fields in the form of "X:"
                section_param_separator = get_next_byte()
                if chr(section_param_separator) != ":" : raise Exception(section_header_errmsg+"missing "":"" separator after parameter")

                # We do not expect compr base values in case of time jump-only headers, so
                # after headers having only one "T:" field.
                if section_param_chr != 'T': expect_compr_base_values = True

                if section_param_chr == "S":
                    station_name = read_string(MAX_HEADER_PARAMETER_NAMES_LEN)
                    # ToDo: check is station name is changed
                elif section_param_chr == "I":
                    instrument_name = read_string(MAX_HEADER_PARAMETER_NAMES_LEN)
                elif section_param_chr == "P":
                    DAQ_app_name = read_string(MAX_HEADER_PARAMETER_NAMES_LEN)
                elif section_param_chr == "D":
                    # Time step definition
                    time_step_str = read_string()
                    if '.' in time_step_str:
                        # Higher timing resolution than 1 sec
                        # E.g.: '0.6' or '1.2'
                        a = time_step_str.index('.')
                        default_time_step__sec =  int(time_step_str [ :a ])
                        subsec_str = time_step_str [ a+1: ]
                        assert len(subsec_str) == 3 # We are prepared only for msec digits
                        default_time_step__msec = int(subsec_str)
                        new_msec_in_first_datacolumn = True
                    else:
                        default_time_step__sec = int(time_step_str)
                        default_time_step__msec = 0
                        new_msec_in_first_datacolumn =  (default_time_step__msec == 0)
                    if msec_in_first_datacolumn is None:
                        msec_in_first_datacolumn = new_msec_in_first_datacolumn
                    else:
                        if msec_in_first_datacolumn != new_msec_in_first_datacolumn: raise Exception("Time-label encoding method changed inside the file!")
                    if (default_time_step__sec<0) or (default_time_step__sec>86400):
                        raise Exception(f"Invalid second step value: {default_time_step__sec}")
                    # time_step__timedelta = datetime.timedelta(days=0, seconds=default_time_step__sec, milliseconds = default_time_step__msec)
                elif section_param_chr == "N":
                    # Nr. of components
                    nr_of_comps_str = read_string()
                    new_nr_of_comps = int(nr_of_comps_str)
                    if nr_of_all_comps is None:
                        # Creating components and reserving buffers:
                        nr_of_all_comps = new_nr_of_comps
                        all_comps = [ ]
                        for i in range(nr_of_all_comps): all_comps = all_comps + [ component_class() ]

                        assert not (msec_in_first_datacolumn is None) # A "D" section should preceed "N" section, so we should know th value of 'msec_in_first_datacolumn'.
                        # If 'msec_in_first_datacolumn' is True, then "all_comps[0]" is the "msec" component
                        stored_comp_indeces = []

                        if comp_indeces is None:
                            # Storing all components
                            for icomp in range(nr_of_all_comps):
                                if (not msec_in_first_datacolumn) or add_msec_to_data or (icomp > 0):
                                    stored_comp_indeces += [icomp]
                        else:
                            # Storing components listed in 'comp_indeces'
                            assert isinstance(comp_indeces, (list,tuple))
                            if msec_in_first_datacolumn and add_msec_to_data:stored_comp_indeces += [0]
                            for comp_index in comp_indeces:
                                comp_index_in_file_comps = comp_index
                                # Shifting component indecex in 'comp_indeces' by one, because file's first component is msec
                                if msec_in_first_datacolumn : comp_index_in_file_comps += 1
                                if comp_index_in_file_comps>=nr_of_all_comps:
                                    raise Exception(f"Component with index {comp_index} is requested, but file has only {nr_of_all_comps} components")
                                stored_comp_indeces += [comp_index_in_file_comps]

                        records = np.empty( [RECORD_BUFFER_INCREASE_STEP, len(stored_comp_indeces)] )
                        assert act_record_index == -1
                        timelabels = np.empty( [RECORD_BUFFER_INCREASE_STEP,1 ] )

                    else:
                        if nr_of_all_comps != new_nr_of_comps: raise Exception("Components number changed within the file ({} -> {})".format(nr_of_all_comps, new_nr_of_comps))

                elif section_param_chr == "C":
                    # Component names
                    if nr_of_all_comps is None: raise Exception("Components number is not defined")
                    comp_names = read_string()
                    comp_names = comp_names.split(',')
                    if len(comp_names) != nr_of_all_comps: raise Exception("Wrong number of component names is defined")
                    for i,comp in enumerate(all_comps):
                        if comp.name == "":
                            comp.name = comp_names[i]
                        else:
                            if comp.name != comp_names[i]:
                                warnings += ["Component nr. {} changed from ""{}"" to ""{}""".format( i+1, comp.name, comp_names[i])]
                                comp.name = comp_names[i]

                    assert not (msec_in_first_datacolumn is None) # A "D" section should preceed "N" section.
                    # If 'msec_in_first_datacolumn' is True, then "all_comps[0]" is the "msec" component
                    if msec_in_first_datacolumn:
                        assert all_comps[0].name == "msec"
                        if (default_time_step__sec == 0) and (all_comps[0].number_total_width == 0):
                            all_comps[0].number_total_width = 3 # all_comps[0] is the millisecond part of time labesl, it can be in the range of 0...999

                elif section_param_chr == "U":
                    # Component units
                    if nr_of_all_comps is None: raise Exception("Components number is not defined")
                    comp_units = read_string()
                    comp_units = comp_units.split(',')
                    if len(comp_units) != nr_of_all_comps: raise Exception("Wrong number of component names is defined")
                    for i,comp in enumerate(all_comps):
                        if comp.unit == "":
                            comp.unit = comp_units[i]
                        else:
                            if comp.unit != comp_units[i]:
                                warnings += ["Component nr. {} unit changed from ""{}"" to ""{}""".format( i+1, comp.unit, comp_units[i])]
                                comp.unit = comp_units[i]

                elif section_param_chr == "T":
                    # if (nr_of_processed_sections>1)  prevTime=mktime(&t);
                    # datetime_str stores time information in the form of %YYYY%MM%DD%hh%mm%ss
                    datetime_str = read_string();
                    record_time = datetime.datetime.strptime(datetime_str, "%Y%m%d%H%M%S")
                    assert record_time.microsecond == 0
                    do_second_step_due_to_msec_overflow = False
                    default_next_msec = 0
                else:
                    raise Exception(f"Unparsed section parameter: {section_param_chr}")

            return expect_compr_base_values

        def add_new_record():
            nonlocal act_record_index, timelabels, records, stored_comp_indeces, record_time, msec_in_first_datacolumn, record_time_msec
            act_record_index = act_record_index + 1
            if len(timelabels)<=act_record_index:
                # Increasing reserved buffers sizes:
                timelabels = np.vstack([timelabels, np.empty([RECORD_BUFFER_INCREASE_STEP, 1]) ])
                records    = np.vstack([records,    np.empty([RECORD_BUFFER_INCREASE_STEP, len(stored_comp_indeces)]) ])
            timelabel = (record_time - EXCEL_DATETIME_EPOCH).total_seconds()/86400.0
            if msec_in_first_datacolumn and (record_time_msec>=0):
                timelabel = timelabel + record_time_msec/86400.0/1000.0
            timelabels[act_record_index] = timelabel

        def process_new_compressed_component_value(negative_sign_flag, delta):
          nonlocal all_comps, comp_index, msec_in_first_datacolumn, default_next_msec, record_time_msec, do_second_step_due_to_msec_overflow

          if negative_sign_flag == 0:
              all_comps[comp_index].compr_base += delta
          else:
              all_comps[comp_index].compr_base -= delta
          decompressed_value = all_comps[comp_index].compr_base / all_comps[comp_index].mult_factor_for_decimals
          if comp_index==0:
              if msec_in_first_datacolumn:
                assert all_comps[comp_index].mult_factor_for_decimals == 1
                record_time_msec =  all_comps[comp_index].compr_base + default_next_msec  # Correction value of 'default_next_msec'
                if (record_time_msec>1000) or (record_time_msec<0):
                    # WARNING(false,"msec value out of range: %d",record_time_msec);
                    # Invalidating msec:
                    record_time_msec=-1
                    default_next_msec=0
                else:
                    if default_time_step__msec>0:
                        # Preparing for the new record:
                        default_next_msec = record_time_msec+default_time_step__msec
                        if default_next_msec >= 1000:
                            assert default_next_msec < 2000
                            do_second_step_due_to_msec_overflow = True
                            default_next_msec=default_next_msec - 1000
                    else:
                        default_next_msec = 0
                        do_second_step_due_to_msec_overflow=False
                decompressed_value = record_time_msec
              add_new_record()

          if comp_index in stored_comp_indeces:
              records[act_record_index, stored_comp_indeces.index(comp_index)] = decompressed_value


        def read_next_component_data():

            nonlocal  nr_of_section_data_records, comp_index, all_comps, nr_of_all_comps, default_time_step__sec, record_time, do_second_step_due_to_msec_overflow
            nonlocal nr_of__0x20_second_time_jump

            b = get_next_byte()

            if (b==0): return False  # new section starting

            if ( (b & 0b11100000) == 0b00100000 ) and ( comp_index==nr_of_all_comps-1):
        	   # Can arrive here before the first component data of each record.
               # Information content: use of "s ssss" second step instead of default_time_step__sec  (for the case of 0 <= sssss < 32)
               # If (default_secstep==0), then this byte is also used to increase the second in subsecond files.

               sec_jump = b & 0b00011111
               nr_of__0x20_second_time_jump += 1
               record_time += datetime.timedelta( seconds = sec_jump )
               do_second_step_due_to_msec_overflow = False
               comp_index = -1
               return True

            comp_index += 1
            if comp_index >= nr_of_all_comps:
                # Increasing time label
                if (default_time_step__sec!=0):
                    record_time += datetime.timedelta( seconds = default_time_step__sec )
                # otherwise the first data columns contains msec value, and second jumps encoded on 0x20 parameter are used
                if do_second_step_due_to_msec_overflow:
                    record_time += datetime.timedelta( seconds = 1 )
                    do_second_step_due_to_msec_overflow = False

                comp_index=0

            if comp_index==0:
                nr_of_section_data_records += 1

            need_decomporess_this_component = (comp_index in stored_comp_indeces) or ( msec_in_first_datacolumn and (comp_index==0) )

            if b==0b00001111:

                 # Invalid component data
                 if comp_index==0:
                     # Can be an invalid time stamp data

                    if (default_time_step__sec==0) and (not msec_in_first_datacolumn):
                		# This component is the abs time - very old format!
                        assert False, "Not supported old format"
                        # t.tm_year=80;t.tm_mon=0;t.tm_mday=1;
                        # t.tm_hour=0;t.tm_min=0;t.tm_sec=0;
                    elif (default_time_step__sec==0) and msec_in_first_datacolumn:
                        # In newer format: an invalid msec data  - but this cannot be happen!
                        # record_time_msec=-1
                        # default_next_msec=0
                        assert False, "invalid millisecond value in file"
                        do_second_step_due_to_msec_overflow=False
                    # Starting new record:
                    add_new_record()

                 # Storing invalid component data
                 if comp_index in stored_comp_indeces:
                     records[act_record_index, stored_comp_indeces.index(comp_index)] = np.nan

                 return True

            if ( (b & 0b10000000)==0b10000000 ):
                # 7 bit delta value
                if need_decomporess_this_component:
                    process_new_compressed_component_value(b & 0x40, b & 0x3F)
                all_comps[comp_index].stat_count[0]+=1
                return True


            if ( (b & 0b01000000)==0b01000000 ):
                b1 = get_next_byte()
                if need_decomporess_this_component:
                    process_new_compressed_component_value(b & 0x20, (b & 0x1F)*0x100 + b1 )
                all_comps[comp_index].stat_count[1]+=1
                return True

            if ( (b & 0b00100000)==0b00100000 ):
        	   # from version vSC2.1
               b1 = get_next_byte()
               b2 = get_next_byte()
               if need_decomporess_this_component:
                   process_new_compressed_component_value(b & 0x10,  (b & 0x0F)*0x10000 + b1*0x100 + b2 )
               all_comps[comp_index].stat_count[2]+=1
               return True

            # From v SC2.1:
            if ( (b & 0b11111100)==0b00001000 ):
               # "0000 10xx" +3 byte formatum
               b1 = get_next_byte()
               b2 = get_next_byte()
               b3 = get_next_byte()
               if need_decomporess_this_component:
                   process_new_compressed_component_value(b & 0x02, (b & 0x01)*0x1000000 + b1*0x10000 + b2*0x100 + b3 )
               all_comps[comp_index].stat_count[3]+=1
               return True

            # "0001 dddd": decimal digits and 4 bytes base value:
            all_comps[comp_index].nr_of_decimals=b & 0x0F
            all_comps[comp_index].mult_factor_for_decimals=10 ** all_comps[comp_index].nr_of_decimals
            b1 = get_next_byte()
            b2 = get_next_byte()
            b3 = get_next_byte()
            b4 = get_next_byte()

            if need_decomporess_this_component:
                l = b4*0x1000000 + b3* 0x10000 +  b2*0x100 +   b1
                if (l & 0x80000000):
                    l = (l ^ 0x80000000) - 0x80000000
                all_comps[comp_index].compr_base = l

                process_new_compressed_component_value(0,0);
            all_comps[comp_index].stat_count[4] += 1
            return True


        file_size = os.stat(file_name)

        if file_size.st_size < 17: raise Exception("File is too small")

        f=open(file_name,"rb")

        file_header = get_next_byte()
        if (file_header != 0): raise Exception("File header error #1")

        EOF = False
        progress_division = 30
        prev_printed_progress_division_index = -1

        while not EOF:
            progress_division_index = int(f.tell() / file_size.st_size * progress_division + 0.5)
            if (progress_division_index != prev_printed_progress_division_index):
                prev_printed_progress_division_index= progress_division_index
                print("-", end='')

            # progress_printer.newline("{:.2f} %".format(relative_filepos))

            expect_compr_base_values  = read_section_params()

            assert not nr_of_all_comps is None
            assert nr_of_all_comps>0


            if expect_compr_base_values:
                for comp in all_comps:
                    comp.compr_base = None

            comp_index = -1
            nr_of_section_data_records=0

            while read_next_component_data():
                if ( index_within_chunk == chunk_len):
                    if ( chunk_len <  CHUNK_SIZE ) or ( f.tell() == file_size.st_size ):
                        EOF = True
                        break

            nr_of_total_data_records += nr_of_section_data_records



    except EOFError:
        errmsg="EOF reached before the file is fully processed"

    except Exception as e:
        errmsg = e.__str__()

    else:
        errmsg = None

    finally:
        # progress_printer.close()
        print("")
        if not f is None: f.close()
        if errmsg is not None: print("Error: "+errmsg)
        if act_record_index>=0:
            # Cutting out the ending unused records
            timelabels = timelabels[:act_record_index+1]
            records    = records[:act_record_index+1,:]
            stored_comps = []
            for comp_index in stored_comp_indeces:
                stored_comps = stored_comps + [ all_comps[comp_index] ]

        else:
            timelabels =None
            records = None
            stored_comps  = None

        return errmsg, warnings, timelabels, records, stored_comps


# Converts an SC file into text file. It does a similar job than 'lscconv' program, with the difference
# that time label here is not written as "hh mm ss msec" format, instead in a numerical format.
def sc_to_ascii_dat(filename):

    errmsg, warnings, timelabels, records, all_comps = sc_load(filename)

    if errmsg is not None:
        print(errmsg)
    else:
        d = np.hstack([timelabels, records])
        delimiter = " "
        comment_sign_before_header = ""

        # First column is the timestamp
        header_line="timestamp"+delimiter

        format_masks = ["%12.6f"]
        for icomp,comp in enumerate(all_comps):

            if (icomp==0) and (comp.name=="msec"):
                header_line += "msec" + delimiter
                format_masks += [ "%03d" ]
            else:
                header_line += comp.name
                if comp.unit !="": header_line +=  "["+comp.unit+"]"
                header_line += delimiter

                format_masks += [ f"%.{comp.nr_of_decimals}f"  ]

        header_line = header_line[:-len(delimiter)]

        np.savetxt("d:\\test.sst", d, fmt = format_masks, header=header_line, delimiter = delimiter, comments = comment_sign_before_header)


# sc_to_ascii_dat("d:\\test.ssc")
