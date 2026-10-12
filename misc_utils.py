# -*- coding: utf-8 -*-
"""

@author: Merényi László

The floating point number encoding date/time is called here as "Excel date".
It is the "OLE Automation date", and it is used also in Delphi/FreePascal and some other software.


"""

import datetime
import numpy as np
import glob
import os
from dateutil.relativedelta import relativedelta



EXCELTIME_ONE_DAY = 1.0
SEC_TO_TIMELABEL_MULT = 1/86400.0  # To converting seconds to Excel date-time value
TIMELABEL_DECIMAL_DIGIT_FOR_SECOND = 5   # Five decimal digits is required for about one second time labels resolution

DEFAULT_DATETIMESTR_FORMATS = [ "%Y.%M.%D %H:%M:%S", "%Y.%M.%D %H:%M"  ]  # Default date/time string parser mask. In amdp project, this is overwritten by 'geomag_utils_preferences.py'


# A 'struct' type variable can be created from this class.
# 'foo = mu.struct()'
# https://stackoverflow.com/questions/35988/c-like-structures-in-python
# 'foo.a = 1'
# ToDo: use NamedTuple or dataclass
#class struct:
#       pass




# Convertes Python time.time() returned datetime struct value into  "Excel" or "Delphi" (or "freepascal") type numerical datetime value.
# This latter one can be handled more easily in Grapher and Excel programs.
# https://stackoverflow.com/questions/9574793/how-to-convert-a-python-datetime-datetime-to-excel-serial-date-number

EXCEL_DATETIME_EPOCH = datetime.datetime(1899, 12, 30)    # Note: not 31st Dec but 30th!
EXCEL_DATETIME_EPOCH_64 = np.datetime64('1899-12-30')
matplotlib_epoch = datetime.datetime(1970, 1, 1)
EXCELTIME_TO_DATETIME_EPOCH_AS_DATETIME = datetime.datetime(1899, 12, 30)  # This gives '0' in LibreOffice calc
EXCELTIME_TO_DATETIME_EPOCH = EXCELTIME_TO_DATETIME_EPOCH_AS_DATETIME.toordinal()


def datetime_2_exceltime(t):
    
    delta = t - EXCEL_DATETIME_EPOCH
    
    if type(t) == np.ndarray:
        delta_td64 = delta.astype('timedelta64[us]')  # or [ns], [s], depending on precision you need
        days_float = delta_td64 / np.timedelta64(1, 'D')
        return days_float
    else:
        return delta.days + delta.seconds / 86400

    # Or with openpyxl: 
    # from openpyxl.utils.datetime import to_excel
    # result = = to_excel(dt)
     
    # if not isinstance(t, (list, np.ndarray)):
    #     assert isinstance(t, datetime.datetime)
    #     delta = t - EXCEL_DATETIME_EPOCH
    #     # Optional: result = delta.total_seconds() / 86400
    #     result = delta.days + delta.seconds / 86400
    #     return 
    # else:
    #     # Converting an array
    #     result = np.empty( len(t) )
    #     for i,_dt in enumerate(t):
    #         result[i] = datetime_2_exceltime(_dt)
    #     return result

def round_to_n_sec(dt, n):
    if isinstance(dt, datetime.datetime):
        # Single conversion
        epoch = datetime.datetime(1970, 1, 1)
        seconds = (dt - epoch).total_seconds()
        rounded = round(seconds / n) * n
        return epoch + datetime.timedelta(seconds=rounded)
    else:
        dt64 = np.asarray(dt, dtype='datetime64[us]')  # elfogad tömböt VAGY egyetlen dátumot is
        epoch = np.datetime64('1970-01-01T00:00:00', 'us')
        
        seconds = (dt64 - epoch) / np.timedelta64(1, 's')   # float másodperc(ek), elemenként
        rounded = np.round(seconds / n) * n                  # np.round() működik tömbön is
        
        rounded_us = np.round(rounded * 1_000_000).astype('int64')
        result = epoch + rounded_us.astype('timedelta64[us]')        
        result = result.tolist()  # datetime64[us] -> datetime conversion
        return result

def exceltime_2_datetime(t, round_to_sec = 1):
    # Other possibiltiy for the conversion:
    # from openpyxl.utils.datetime import from_excel
    # dt = from_excel(46245.605636574)
    
    if type(t) == np.ndarray:
        if type(round_to_sec)==int:
          # Conversion with one second resolution:
          seconds = np.round(t * 86400).astype('int64')
          deltas = seconds.astype('timedelta64[s]')
          result = EXCEL_DATETIME_EPOCH_64 + deltas    
          result = result.tolist()  # datetime64[us] -> datetime conversion
          if round_to_sec != 1:
               result = round_to_n_sec(result, round_to_sec)
        else:
          # Conversion with sub-second resolution
          assert round_to_sec is None
          micros = np.round(t * 86400 * 1_000_000).astype('int64') # Napok -> mikroszekundum pontosságú egész számok, majd timedelta64-gyé alakítva
          deltas = micros.astype('timedelta64[us]')
          result = EXCEL_DATETIME_EPOCH_64 + deltas
          result = result.tolist()  # datetime64[us] -> datetime conversion
    else:
        result = EXCEL_DATETIME_EPOCH + datetime.timedelta(days=t)
        # days, remaining_time = divmod(t, 1)
        # seconds, microseconds = divmod( remaining_time * 86400.0 , 1)
        # microseconds = microseconds * 1e6
        if type(round_to_sec)==int:
            #if microseconds>500000:
            #    seconds += 1
            #microseconds = 0
            if result.microsecond >= 500_000:
                result = result + datetime.timedelta(seconds=1)
            result = result.replace(microsecond=0)
            if round_to_sec != 1:
                result = round_to_n_sec(result, round_to_sec)
        else:
            assert(round_to_sec is None), f"{round_to_sec}"
    
        # timedelta returns elapsed time since an EPOCH
        # delta_from_epoch = datetime.timedelta( days, seconds, microseconds )
    

    return result

# ToDo: unify with 'exceltime_2_datetime'
# def exceltime_to_datetime(excel_date, round_to = "sec"):
#     # From 'https://stackoverflow.com/questions/31359150/convert-date-from-excel-in-number-format-to-date-format-python'
#     # In that solution, EPOCH is 1900.01.01, but '2' is removed explicitly
#     def floatHourToTime(fh):
#         h, remainder = divmod(fh*24, 1)
#         m, remainder = divmod(remainder*60, 1)
#         sec, remainder = divmod(remainder*60, 1)
#         msec = remainder  *1000
#         if msec > 999: msec = 999
#         return int(h),int(m),int(sec),int(msec)

#     if round_to=="sec":
#         excel_date = np.round  (excel_date * 86400.0) / 86400.0
#     elif round_to=="min":
#         excel_date = np.round  (excel_date * 1440.0)  / 1440.0
#     elif round_to=="msec":
#         excel_date = np.round  (excel_date * 86400.0 * 1000.0) / (86400.0 *1000.0)
#     else:
#         assert False

#     if isinstance(excel_date, np.ndarray):
#         assert False, "list/ndarray is not supported here -> use 'exceltime_2_datetime()' instead"
#         dt = datetime.datetime.fromordinal(EXCELTIME_TO_DATETIME_EPOCH + excel_date.astype(int) )
#     else:
#         dt = datetime.datetime.fromordinal(EXCELTIME_TO_DATETIME_EPOCH + int(excel_date) )

#     hour, minute, second, msec = floatHourToTime(excel_date % 1)

#     if round_to=="sec":
#         dt = dt.replace(hour=hour, minute=minute, second=second, microsecond = 0)
#         if msec > 500: dt = dt + datetime.timedelta(seconds=1)
#     elif round_to=="min":
#         dt = dt.replace(hour=hour, minute=minute, second=0, microsecond = 0)
#         if second >= 30:dt = dt + datetime.timedelta(minutes=1)
#     elif round_to=="msec":
#         microsec = msec * 1000.0
#         assert round(microsec)<1e6
#         dt = dt.replace(hour=hour, minute=minute, second=0, microsecond = round(microsec))
#     else:
#         assert False
#     # timetuple: tt = dt.timetuple()

#     if (round_to == "sec"):
#         check = exceltime_2_datetime(excel_date, "sec")
#         diff = dt-check
#         if isinstance(excel_date, float) : diff = [ diff ]
#         for idiff in diff:
#             assert abs( idiff.total_seconds()  )<1.0
#     return dt


def datetimestr_2_datetime(datetime_str, format = None):
     if format is None: 
         # Using default date-time formats
         result = None         
         for format in DEFAULT_DATETIMESTR_FORMATS: 
             try:
                 result = datetime.datetime.strptime(datetime_str,format)
                 break  # Using first format that fits to 'datetime_str'
             except ValueError:
                 pass
         if result is None: 
              raise ValueError("No default mask to parse  \"{}\" date-time or wrong date-time".format(datetime_str))
     else:
         result = datetime.datetime.strptime(datetime_str,format)
     return result
    

def datetimestr_2_exceltime(datetime_str, format = None):
    result = datetimestr_2_datetime(datetime_str, format)
    if result is None: return result
    return datetime_2_exceltime(result)

def exceltime_2_datetimestr(excel_time, format = ""):  # Format can be: "%Y.%m.%d %H:%M:%S"
     datetime=exceltime_2_datetime(excel_time)
     if format == "":
         # return "{}".format(datetime)  # Using Python default datetime-to-string format
         return datetime.strftime(DEFAULT_DATETIMESTR_FORMATS[0])
     else:
         return datetime.strftime(format)


def exceltime_2_matplotlib_time(d):
    epoch_correction = EXCELTIME_TO_DATETIME_EPOCH_AS_DATETIME - matplotlib_epoch
    if isinstance(d,list):
        return [ i + float(epoch_correction.days) for i in d]
    return  d + float(epoch_correction.days)

def matplotlib_time_2_exceltime(d):
    epoch_correction = EXCELTIME_TO_DATETIME_EPOCH_AS_DATETIME - matplotlib_epoch
    if isinstance(d, list):
        return [i - float(epoch_correction.days) for i in d]
    return d - float(epoch_correction.days)


def change_percent_get(current, previous):
    if current == previous:
        return 0
    try:
        return (abs(current - previous) / previous) * 100.0
    except ZeroDivisionError:
        return float('inf')

# From https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
def rotate(around_axis, angle, vec):
    res = [0,0,0]
    if around_axis=='X':
        res[0] = vec[0]
        res[1] = vec[1]*np.cos(angle) - vec[2]*np.sin(angle)
        res[2] = vec[1]*np.sin(angle) + vec[2]*np.cos(angle)
    elif around_axis=='Y':
        res[0] = vec[0] * np.cos(angle) + vec[2]*np.sin(angle)
        res[1] = vec[1]
        res[2] =-vec[0] * np.sin(angle) + vec[2]*np.cos(angle)
    elif (around_axis=='Z'):
        res[0] = vec[0] * np.cos(angle) - vec[1]*np.sin(angle)
        res[1] = vec[0] * np.sin(angle) + vec[1]*np.cos(angle)
        res[2] = vec[2]
    else:
        assert(False)
    return res

def threeaxis_rotate(x_angle, y_angle, z_angle, vec):
    # X axis rot:
    xres = [vec[0], vec[1]*np.cos(x_angle) - vec[2]*np.sin(x_angle), vec[1]*np.sin(x_angle) + vec[2]*np.cos(x_angle) ]
    # Y axis rot:
    xyres = [ xres[0] * np.cos(y_angle) + xres[2]*np.sin(y_angle), xres[1], -xres[0] * np.sin(y_angle) + xres[2]*np.cos(y_angle) ]
    # Z axis rot:
    xyzres = [ xyres[0] * np.cos(z_angle) - xyres[1]*np.sin(z_angle), xyres[0] * np.sin(z_angle) + xyres[1]*np.cos(z_angle), xyres[2] ]
    return xyzres



def find_outliers(da, limit, datetime):
    outlier_indeces = []
    da_copy = np.copy(da)
    while True:
        max_index = np.argmax(da_copy)
        if (da_copy[max_index] < limit): break  # No more outlier
        print('Outlier: index={:d} value={:.3f}  time={:.6f} '.format(max_index, da[max_index], datetime[max_index]) )
        outlier_indeces = outlier_indeces  + [max_index]
        da_copy[max_index] = 0
    return outlier_indeces


# 'mask' formats: "string", "*string", "string*"
# If 'mask' is a regular string: returns 'mask == s'
# If 'mask' starts with an asterisk, then returns 'true' if 's' ends with the remining part of the mask
# If 'mask' ends with an asterisk, then returns 'true' if 's' starts with the starting part of the mask
def compare_with_string_mask(mask, s):
    if isinstance(mask, bytes): mask = mask.decode()
    if isinstance(s, bytes): s = s.decode()

    mask_len = len(mask)
    if mask_len > 0:
        if mask[0] == '*':
            if mask == "*": return True
            return s[ -(mask_len-1) : ] == mask[1:]
        if mask[mask_len-1] == '*':
            return s[ : mask_len-1 ] == mask[ : mask_len-1]
    return mask == s

# Linear fit
def LSQ(x,y, handle_nans = False, return_resid_stddev = False, return_resid_array = False):

    missing = np.isnan(x)
    if np.any(missing):
        assert handle_nans,"There are NAN's in x data"
        x=x[~missing]
        y=y[~missing]
    missing = np.isnan(y)
    if np.any(missing):
        assert handle_nans,"There are NAN's in y data"
        x=x[~missing]
        y=y[~missing]

    assert len(x)>2
    assert len(x) == len(y)

    A = np.vstack( [x, np.ones(len(x))] ).T
    linreg_m, linreg_c = np.linalg.lstsq(A, y, rcond=None)[0]
    
    resid_info = []
    if return_resid_array or return_resid_stddev:
        resid = y - (linreg_c + linreg_m * x)
        if return_resid_stddev:
            resid_info += [ np.std(resid) ] 
        if return_resid_array:
            resid_info += [ resid ] 
    
    return  linreg_c, linreg_m, resid_info

# If 'p' is a simple item ( not a list and not an np.darray), then this function
# will create and return a list with a simple element of 'p'.
# Otherwise, 'p' is returned withouth conversion
def ensure_list(p):
    if isinstance(p, (list, np.ndarray, tuple)): return p
    # if p is None: return [] => if None should be an empty list
    return [ p ]

def ensure_file_name_extension(filename,ext):
    if ext[0] != ".": ext = "."+ext
    ext_len=len(ext)
    if (filename[-ext_len:] != ext):  filename = filename + ext
    return filename


# Formats a string representation of a list/tuple
# Used to output lists to user.
def list_to_string(l, separator=","):
    assert isinstance(l, (tuple,list))
    result = ""
    for i,item in enumerate(l):
        if i>0: result = result + separator
        if isinstance(item,bytes):
            result = result + item.decode('UTF-8')
        elif isinstance(item,int):
            result = result + str(item)
        elif isinstance(item,str):
            result = result + item
        else:
            assert False, "Not yet implemented"

    return result

def plural_suffix(list_or_tuple_or_array_or_number, singular_result = "", plural_result = "s"):
    if isinstance( list_or_tuple_or_array_or_number, ( int, np.int64) ): 
        count = list_or_tuple_or_array_or_number
    else:
        count = len(list_or_tuple_or_array_or_number)
    return plural_result if count > 1 else singular_result


# Checks if there is at least one NAN in an entire structued array
def structured_array_has_nans(structed_array):
    for comp in structed_array.dtype.names:
        if np.any( np.isnan(structed_array[comp])):
            return True
    return False
    # Alternatice way: has_nan = any(np.isnan(a[name]).any() for name in a.dtype.names  if a.dtype[name].kind == 'f' )

def list_timelabels(msg, timelabels, max_printed_lines = 10):
    print(msg)
    for i in range(min(max_printed_lines,len(timelabels) ) ):
        t = timelabels[i]
        print("\t{:.6f} {:}".format(t, exceltime_2_datetime(t)) )
    if len(timelabels) > max_printed_lines:
        print("... (there are {} more timelabels))".format(len(timelabels) - max_printed_lines))


# Tries to determine a format mask for an array
def auto_format_mask(r, comp_name, mask_type = "python-default"):

    if (isinstance(r, (np.str_))) or (r.dtype.kind in {'U', 'S'}):
        return "%s"

    if np.all( np.isnan(r) ):
        # Component has only nan values, that will be written as "nan"
        return "%f"
        
    minval = np.nanmin(r)
    maxval = np.nanmax(r)
    abs_max_val = max(abs(minval),abs(maxval))

    diff = maxval-minval

    if diff>=10.0:
        digits = 2
    elif diff>=1.0:
        digits = 6
    elif diff>=0.1:
        digits = 7
    else:
        digits = 6

    if mask_type == "python-default":
        return "{{:.{}f}}".format(digits)
    if mask_type == "numpy-fix-len":
        # Python 'str(xxx)' probably auto determin decimal number required for correct representation
        s = str(abs_max_val)
        if not '.' in s:
            # Integer number 
            int_part_len = len(s)
        else:
            int_part_len = s.index('.')
        total_field_len = int_part_len + 1 + digits
        if minval<0: total_field_len += 1
        return "%{}.{}f".format(total_field_len, digits)

    assert False, "Unsupported mask_type: {}".format(mask_type)



# Extension of numpy.savetxt() with some extra options. 
# It saves a numpy structured array into a text file, with header and auto
# calcualted format masks.
#  'compoment_parameters' can encode parameters in this formats:
# 1. 'compoment_parameters' is a dictinoary, with component names as keys.
#       Values are struckt-like class instances, with a 'format_mask' and 'unit' member variables.
def savetxt(filename, r, component_parameters = None):
    
    assert isinstance(r, np.ndarray)
    assert r.ndim == 1
    
    delimiter = " "
    comment_sign_before_header = ""

    format_masks = []
    header_line = ""

    for icomp,compname in enumerate(r.dtype.names):
        
        if (component_parameters is None) or (component_parameters[compname].unit == ""):
            component_name_in_header = compname
        else:
            component_name_in_header = compname+'['+component_parameters[compname].unit+']'
        
        header_line += component_name_in_header + delimiter
        
        if component_parameters is not None:
            format_masks += [ component_parameters[compname].format_mask  ]         
        elif r.dtype.kind in {'U', 'S'}:
            format_masks += [ "%s" ] # string field
        elif r.dtype[icomp].name[:3] =='int':
            if compname in ['hour','min','sec']:
                format_masks += [ '%02d' ]
            else:
                format_masks += [ '%d' ]
        else:
            format_masks += [ auto_format_mask(r[compname], compname, "numpy-fix-len") ]

    header_line = header_line[:-len(delimiter)]

    np.savetxt(filename, r, fmt = format_masks, header=header_line, delimiter = delimiter, comments = comment_sign_before_header)


# Returns a list of days, between two dates (start and end dates are included)
# From: https://stackoverflow.com/questions/1060279/iterating-through-a-range-of-dates-in-python
def day_range(start_date: datetime.date, end_date: datetime.date):
    days = int((end_date - start_date).days) + 1
    for n in range(days):
        yield start_date + datetime.timedelta(n)



def datetime_step(step,time_unit):
    if time_unit == "min":
        return step/1440.0
    assert False, "not implemented"


# If the parameter has a date/time mask, then this function extract it to real, existing file names,
# that fit to the mask. 
# 'filenamsmask' parameter can be a single string (with optionally mask elements) or a list/tuple of filenames/masks.
# 'base_path' is a static text, without mask.
# Returned value is a list of fitted filenames, and a list of used file-masks
def filenamemask_to_filenames(base_path, filenamsmask, datetime_filters = None):
    
    if isinstance(filenamsmask, (list, tuple)): 
        # List of masks (e.g. ADS raw files definition)
        all_filenames = [] 
        all_file_masks = []
        for act_mask in filenamsmask: 
            new_filenames, new_filemask = filenamemask_to_filenames(base_path, act_mask, datetime_filters)
            assert isinstance(new_filenames, (list,tuple))
            all_filenames += new_filenames 
            assert len(new_filemask) == 1 # It should be a list with one element
            if not new_filemask[0] in all_file_masks: all_file_masks += new_filemask
        return all_filenames, all_file_masks
    
    
    #if not '%' in filenamsmask: 
    #    return [ [filenamsmask], [filenamsmask] ] 
    #    # 'mask' is an explicit file name. 
    #    # Todo: check if file is exists
    

    # Making dir pattern:
    dir_mask = filenamsmask
    file_interval = None
    
    if "%d" in dir_mask:
        file_interval = relativedelta(days = 1)
        dir_mask = dir_mask.replace("%d","??")
    if "%m" in dir_mask:
        if file_interval is None: file_interval = relativedelta(months = 1)
        dir_mask = dir_mask.replace("%m","??")
    if "%Y" in dir_mask:
        if file_interval is None: file_interval = relativedelta(years = 1)
        dir_mask = dir_mask.replace("%Y","????")
        
    # dir_mask = mask.replace("%Y","????").replace("%y","??").replace("%m","??").replace("%d","??")
    result_file_names = []    
    result_file_mask = [ dir_mask ] 
    if base_path is None: base_path = ""  # 'base_path' is a constans string, it can not have wildcards
    dir_mask = os.path.join(base_path,dir_mask)
    filelist = glob.glob(dir_mask)
    for filename in filelist:
        # filename = filename[ len(base_path): ]  do not work always!
        filename = os.path.basename(filename)
        if file_interval is not None:
            # File name encodes a date
            try: 
                date_from_filename = datetime.datetime.strptime(filename, filenamsmask)
                # We get here if there is no exception, so file name fits to the mask
                # if base_path != "":  filename = os.path.basename(filename)
                if datetime_filters is not None:
                    if "FROM-DATE" in datetime_filters.keys():
                        from_date_def = datetime_filters["FROM-DATE"]
                        if isinstance(from_date_def, (list,tuple)):
                            from_date_datetime = from_date_def[0]  # Getting the 'datetime' representation (first element of the tuple)
                        else:
                            from_date_datetime = from_date_def   # It should be datetime
                        if date_from_filename + file_interval <= from_date_datetime: continue
                    if "TO-DATE" in datetime_filters.keys():
                        to_date_def = datetime_filters["TO-DATE"]
                        if isinstance(to_date_def, (list,tuple)):
                            to_date_datetime = from_date_def[0]  # Getting the 'datetime' representation (first element of the tuple)
                        else:
                            to_date_datetime = to_date_def  # It should be datetime
                        if date_from_filename > to_date_datetime: continue
                    
                
                result_file_names += [filename]
            except ValueError:
                pass
        else:
            result_file_names += [filename]
    return result_file_names, result_file_mask


def mask_to_date_interval(mask):
    if "%d" in mask: return "day"
    if "%m" in mask: return "month"
    if "%Y" in mask: return "year"
    assert False
        
def generate_date_range(start_date, end_date, interval_type="day", step=1):
    # Convert strings to datetime objects
    b_use_exceltime = isinstance(start_date, (float,int))
    if b_use_exceltime:
        start_date = exceltime_2_datetime(start_date)
        end_date = exceltime_2_datetime(end_date)
        
    # start_date = datetime.strptime(start_str, "%Y-%m-%d")
    # end_date = datetime.strptime(end_str, "%Y-%m-%d")
    
    start_date = start_date.date()
    end_date  = end_date.date()    
    
    # Configure the step increment
    if interval_type == "day":
        delta = datetime.timedelta(days=step)
    elif interval_type == "month":
        start_date = start_date.replace( day = 1 )
        end_date   = end_date.replace( day = 1 )
        delta = relativedelta(months=step)        
    elif interval_type == "year":
        start_date = start_date.replace( day = 1, month = 1 )
        end_date   = end_date.replace( day = 1, month = 1 )
        delta = relativedelta(years=step)
    else:
        raise ValueError("interval_type must be 'day', 'month', or 'year'")
        
    # Generate the date list
    date_list = []
    current_date = start_date
    
    while current_date <= end_date:
        current_datetime = datetime.datetime.combine( current_date, datetime.time() )        
        if b_use_exceltime:
            date_list.append( datetime_2_exceltime(current_datetime) )
        else:
            date_list.append(current_datetime)
        current_date += delta
        
    return date_list    
    

