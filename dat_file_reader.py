# -*- coding: utf-8 -*-
"""

dat_file_reader.py

Reading time stamped data records from data files.
Supported file formats:
- ASCII DAT files (including IAGA-2002 format files)
- maglin saved SC compressed binary files
- CDF (currently with some hard wired conversions)

Created on Jul 18 2021
@author: Merényi László

"""

import os
import numpy as np
from numpy.lib import recfunctions as rfn
import datetime
import re  # For interval parsing
from enum import Enum

import misc_utils as mu
import sc_file_reader
from collections import OrderedDict


SUPPORTED_FILE_FORMATS =  ( "ASCII-dat", "SC", "CDF", "IAGA-2002")

CUSTOM_TIMESTAMP_MASKS = [ "%A", "%ms" ]

class file_read_result_data_type(Enum):
     records_for_timestamps = 0      # Either raw records, or interpolated records, depending on 
     _1min_arithm_avg = 1
     arithm_avg_within_period = 2
     all_within_period = 53


TIME_COMPONENT_NAME__HOUR      = "__hour__"
TIME_COMPONENT_NAME__MINUTE    = "__min__"
TIME_COMPONENT_NAME__SECOND    = "__sec__"
TIME_COMPONENT_NAME__MILLISECOND = "__msec__"
TIME_COMPONENT_NAME__EXCELTIME = "__A__"
    

# "AUTO_FOR_SECOND" is only supported by 'ads_utils.load_mag_data()': it changes to exact_time_match or to file_read_method_enum.raw_interpolated, depending on srouce file interval (sec or min)

# _NO_DATA = []  # This is returned by some function if no data could be read for a datetime
DEF_MIN_SAMPLES_FOR_AVG_CALC = 3
MAX_BUFFERED_SIZE = 200 * 1_000_000  # To store previous read data files in memory, when a new file is open 200 MByte
warn_zero_timeindeces_diff = False

CACHE_DIR=None 
CACHE_FILE_EXT="npy"


TIME_UNITS = {
    "ms": "milliseconds",
    "msec": "milliseconds",
    "millisec": "milliseconds",
    "milliseconds": "milliseconds",    
    "sec": "seconds",
    "s": "seconds",
    "min": "minutes",
    "m": "minutes",
    "hour": "hours",
    "h": "hours",
    "day": "days",
    "d": "days",
}

def find_closest_index(sorted_arr, value):
    idx = np.searchsorted(sorted_arr, value)
    
    # handle edge cases: value is beyond either end
    if idx == 0: return 0
    if idx == len(sorted_arr): return idx-1
    
    # compare the two neighbors
    before = sorted_arr[idx - 1]
    after = sorted_arr[idx]
    if after - value < value - before:
        return idx
    else:
        return idx-1
    

def parse_interval(s):
    
    if s.endswith("/not-aligned"):
        s = s[:-12]
        aligned = False
    else:
        if s.endswith("/aligned"): s = s[-8]
        aligned = True

        
    match = re.match(r"^(\d+)\s*([a-zA-Z]+)$", s.strip())
    if not match:
        raise ValueError(f"Interval parsing failed: {s!r}")
    time_value, time_unit = match.groups()
    time_unit = time_unit.lower()
    if time_unit not in TIME_UNITS:
        raise ValueError(f"Unsupported interval unit: {time_unit!r}")
    time_unit = TIME_UNITS[time_unit]
    
    datetime_timedelta = datetime.timedelta(**{time_unit: int(time_value)})
    
    if time_unit == "seconds":
        interval_sec = int(time_value)
        assert datetime_timedelta.total_seconds() == interval_sec
    elif time_unit == "minutes" :
        interval_sec = int(time_value) * 60
        assert datetime_timedelta.total_seconds() == interval_sec
    elif time_unit == "hours":
        interval_sec = int(time_value) * 3600
        assert datetime_timedelta.total_seconds() == interval_sec
    else:
        assert False
    
    return datetime_timedelta, interval_sec, time_unit, aligned

class dat_file_reader_class:

    # 'range_filter' must be a dict of 'min' and 'max':
    # 'min'/'max' Items can be:
    #   1. dict of {compname:value}
    #   2. list of values, for each components

    def __init__(self, path, file_name_mask, file_record_interval_def, components, data_comps_file_colindeces = None, conversion_defs = None, FILE_COLINDEX_HOUR=None, file_format="", data_descr = "", range_filter = {'min':None, 'max':None}):

        self.file_format = file_format        
        self.data_descr = data_descr
        
        self.path = os.path.join(path,'')
        self.path_has_datetime_mask = "%" in self.path        
        self.file_name_mask = file_name_mask  # A string with formatting mask for year/month/day
        
        self.FILE_COLINDEX_HOUR = None
        self.FILE_COLINDEX_MIN  = None
        self.FILE_COLINDEX_SEC  = None
        self.FILE_COLINDEX_MSEC = None
        self.FILE_COLINDEX_EXCELDATETIME = None
        
        self.info_lines_prefix = ""
        self.detailed_file_open_line = True
        
        assert isinstance(components, (list,tuple))
        self.comp_names = []
        self.data_comps_file_colindeces = [] if data_comps_file_colindeces is None  else data_comps_file_colindeces  # In 'none' this it will be automatic (first data component is after the last timelabel)                
        # 'components' can hold data components only, or cobined timestamp elements + data components

        auto_determine_data_comp_file_colindeces =  data_comps_file_colindeces is None # If not supplied as parameter: by default, we assume that component order in 'components' does not correspond to component order in file. 
        self.time_components_file_colindeces = []
        self.time_components_names = []
        act_column_index = 0

        auto_determine_timestamp_components = True          
        timestamp_component_names = [ "%Y", "%m", "%d", "%H", "%M", "%S" ] + CUSTOM_TIMESTAMP_MASKS
        for comp in components:
            if comp in timestamp_component_names:
                auto_determine_timestamp_components = False # As timestamp components are explicitly defined, we do not determine it automatically
                break
        
        # 'components' can be just a simple enumuration of the data compoments. In that case file structure is determined either automatically, or from 'data_comps_file_colindeces' parameter
        # 'components' can also be a full record definition: the enumeration of timestamp components, plus eventual skipped columns, plus compoments (used data columns)
        for comp in components:
            if comp[0] != '%':
                # Data component
                self.comp_names += [ comp ]
                if auto_determine_data_comp_file_colindeces and (not auto_determine_timestamp_components):
                    self.data_comps_file_colindeces += [ act_column_index ]
            elif comp == '%skip':
                # if data_comps_file_colindeces is None: auto_determine_data_comp_file_colindeces = True # We assume, that data components follow the column order
                pass
            elif comp in timestamp_component_names:
                # A timestamp component
                if data_comps_file_colindeces is None: auto_determine_data_comp_file_colindeces = True # We assume, that data components follow the column order
                self.time_components_file_colindeces += [ act_column_index ]
                match comp:
                    case "%H":
                        assert self.FILE_COLINDEX_HOUR is None, "Duplicate hour definition"                        
                        self.FILE_COLINDEX_HOUR = act_column_index
                        if FILE_COLINDEX_HOUR is not None: assert FILE_COLINDEX_HOUR == self.FILE_COLINDEX_HOUR, "Mismatch in hour timelabel definition!"
                        self.time_components_names += [ TIME_COMPONENT_NAME__HOUR ]
                    case "%M":
                        assert self.FILE_COLINDEX_MIN is None, "Duplicate minute definition"                        
                        self.FILE_COLINDEX_MIN = act_column_index                
                        self.time_components_names += [ TIME_COMPONENT_NAME__MINUTE ]
                    case "%S":
                        assert self.FILE_COLINDEX_SEC is None, "Duplicate second definition"                        
                        self.FILE_COLINDEX_SEC = act_column_index                
                        self.time_components_names += [ TIME_COMPONENT_NAME__SECOND ]
                    case "%ms":                        
                        assert self.FILE_COLINDEX_MSEC is None, "Duplicate msec definition"                        
                        self.FILE_COLINDEX_MSEC = act_column_index                
                        self.time_components_names += [ TIME_COMPONENT_NAME__MILLISECOND ]
                    case "%A":
                        assert False, "Feature not yet supported"
                        assert self.FILE_COLINDEX_EXCELDATETIME is None, "Duplicate excel date/time definition"                        
                        self.FILE_COLINDEX_EXCELDATETIME = act_column_index                
                        self.time_components_names += [ TIME_COMPONENT_NAME__EXCELTIME ]
                    case _:
                        raise ValueError("Not yet supported timestamp mask: {}".format(comp))
            act_column_index += 1
            
        if len(self.data_comps_file_colindeces) > 0:
            auto_determine_data_comp_file_colindeces = False # No need for using an component colindex auto determation (see later), as it was determined from the component list
            
        self.nof_comps = len(self.comp_names)
        assert self.nof_comps > 0 , "No data column defined"
        
        self.file_record_interval_timedelta, self.file_record_interval_sec, self.file_record_interval_unit, self.file_record_interval_aligned = parse_interval(file_record_interval_def)
        
        if (self.file_record_interval_unit == "seconds") and (self.file_record_interval_aligned): 
            self.timestamp_index_divider_sec = self.file_record_interval_sec
            # If file record interval is more than 1 sec (e.g. GSM magnetometer data), 
            # then take care in case of not exact second positions!
        
        if self.file_record_interval_timedelta.microseconds  !=0:
            self.timeindex_interval_in_sec = self.file_record_interval_timedelta.total_seconds()
        elif self.file_record_interval_timedelta.days == 0:
            self.timeindex_interval_in_sec = self.file_record_interval_timedelta.seconds  # could be named also as filerecord_interval_sec
            assert type(self.timeindex_interval_in_sec) == int
        else:
            assert False, "Day interval is not yet supported"
        
        if "%d" in file_name_mask:
            assert "%m" in file_name_mask, file_name_mask
            assert "%Y" in file_name_mask, file_name_mask
            self.one_file_time_period_sec = 86400.0
            self.file_name_interval_type = "days" # Recently  only daily files are supporte
        elif "%m" in file_name_mask:
            assert "%Y" in file_name_mask
            self.one_file_time_period_sec = None  # Not yet supported, will fail later
            self.file_name_interval_type = "months" 
        elif "%Y" in file_name_mask:
            self.one_file_time_period_sec = None  # Not yet supported, will fail later
            self.file_name_interval_type = "years" 
            
        if auto_determine_timestamp_components:
            self.FILE_COLINDEX_HOUR = 0 if FILE_COLINDEX_HOUR is None else FILE_COLINDEX_HOUR # By default, it is the first column
            self.FILE_COLINDEX_MIN  = self.FILE_COLINDEX_HOUR + 1
            
            match self.file_record_interval_unit:
                case "milliseconds":
                    self.FILE_COLINDEX_SEC  = self.FILE_COLINDEX_MIN + 1
                    self.FILE_COLINDEX_MSEC = self.FILE_COLINDEX_SEC + 1
                    self.time_components_file_colindeces = [ self.FILE_COLINDEX_HOUR, self.FILE_COLINDEX_MIN, self.FILE_COLINDEX_SEC, self.FILE_COLINDEX_MSEC ]
                    self.time_components_names = [ TIME_COMPONENT_NAME__HOUR, TIME_COMPONENT_NAME__MINUTE, TIME_COMPONENT_NAME__SECOND, TIME_COMPONENT_NAME__MILLISECOND ] 
                case "seconds":
                    self.FILE_COLINDEX_SEC = self.FILE_COLINDEX_MIN + 1
                    self.time_components_file_colindeces = [ self.FILE_COLINDEX_HOUR, self.FILE_COLINDEX_MIN, self.FILE_COLINDEX_SEC ]
                    self.time_components_names = [ TIME_COMPONENT_NAME__HOUR, TIME_COMPONENT_NAME__MINUTE, TIME_COMPONENT_NAME__SECOND ] 
                case "minutes":
                    self.time_components_file_colindeces = [ self.FILE_COLINDEX_HOUR, self.FILE_COLINDEX_MIN ]
                    self.time_components_names = [ TIME_COMPONENT_NAME__HOUR, TIME_COMPONENT_NAME__MINUTE] 
                case _:
                    raise ValueError("Unsupported file interval type: {}".format( self.file_record_interval_unit) )

        if auto_determine_data_comp_file_colindeces:
            # Automatic data column indexes
            if type(self.timeindex_interval_in_sec) == float:
                # Sub-second - assuming milliste
                firstcomp_colindex = self.FILE_COLINDEX_HOUR + 4
            else:
                firstcomp_colindex = self.FILE_COLINDEX_HOUR + 2  # assuming hour+minute are always in the fole
                if (self.timeindex_interval_in_sec < 60 ) or  ((self.record_interval_sec % 60  ) != 0):
                    firstcomp_colindex += 1 # Second component
            self.data_comps_file_colindeces = list( range(firstcomp_colindex, firstcomp_colindex + len(self.comp_names)) )
        else:
            assert isinstance(self.data_comps_file_colindeces, list)
            
        assert self.nof_comps == len(self.data_comps_file_colindeces)
        
        self.leave_partial_records = True
        self.fix_header_lines = 1
        self.load_file_name = ""
        self.not_found_file_names = []
        self.reported_missing_files = []
        self.load_file_name_base = ""
        self.interpolated_data_result_halfwindow_size = None
        self.range_filter = range_filter.copy() 
        self.conversion_defs = conversion_defs
        self.no_record_found_return_data=np.empty([1, self.nof_comps]) 
        self.no_record_found_return_data[:] = np.nan
        self.file_record_comments = '#'

        # We buffer the content of the last several read files - useful when running averages / filtering are used that passes day limits
        self.buffered_files = OrderedDict()   # A dict of dict, in this format: { filename: { 'size':nnn, 'data':buff-array (can be an empty array),  'timeindeces':timeindeces, 'timestamps':timelabels}  } 
        self.sum_buffer_data_size = 0
        
        self.set_processing_parameters()
        


    def file_name_get(self, dt_stru):

        if '{:' in self.file_name_mask:
            # Year-month-day order, with Python format mask codes (e.g. {:04d}
            return self.file_name_mask.format(dt_stru.year,dt_stru.month,dt_stru.day)
        else:
            return dt_stru.strftime(self.file_name_mask)


    def calc_timestamp_indeces(self, hour, minute, second, msec=None):
        
        if not self.file_record_interval_aligned:
            return None
        if self.file_record_interval_unit == "milliseconds":
            # ToDo: ensure that interval is 1 msec
            timestamp_indeces = np.round(hour*3600*1000 + minute*60*1000 + second*1000 + msec).astype(int)
        if self.file_record_interval_unit == "seconds":
            assert type(self.timestamp_index_divider_sec) == int
            timestamp_indeces = np.round( (hour*3600 +  minute*60 + second) // self.timestamp_index_divider_sec ).astype(int)
        elif self.file_record_interval_unit == "minutes":
            timestamp_indeces = np.round(hour*60 +  minute).astype(int)
        else:
            assert False
        return timestamp_indeces

    
    def calc_timestamps(self, date_from_filename, hour, minute, second, msec=None):
        if self.file_record_interval_unit == "milliseconds":
            timestamps =  date_from_filename + (hour*3600 +  minute*60 + second + msec/1000.0 ) / 86400.0
        if self.file_record_interval_unit == "seconds":
            timestamps = date_from_filename + (hour*3600 +  minute*60 + second ) / 86400.0
        elif self.file_record_interval_type == "minutes":
            timestamps = date_from_filename + (hour*3600 +  minute*60 ) / 86400.0
        else:
            assert False
        return timestamps
    

    # For daily files, it removes hour/minute/second from a time stamp, to get a time stamp representing the file name
    # Works both for excel datetime values and for datetime.datetime struct
    def  truncate_datetime_to_filename_period(self, timelabel):
        # For date/time truncation hints, see 'get_minute_timelabels()'
        assert self.file_name_interval_type == "days"
        if isinstance(timelabel,datetime.datetime):
            # Excel datetime is received
            return timelabel.replace(hour=0, minute=0, second=0, microsecond=0)
        assert isinstance(timelabel,float)
        return float( int(timelabel) )

    def save_cached_file(self,cached_file_name):
        if CACHE_DIR is None: return
        assert self.load_file_timestamp_indeces is not None, "Fix this"
        self.load_file_timestamp_indeces.shape = ( len(self.load_file_timestamp_indeces), 1  )
        self.load_file_timestamps.shape = ( len(self.load_file_timestamps), 1  )
        stacked_data = np.hstack([self.load_file_timestamp_indeces, self.load_file_timestamps, self.load_file_data])
        np.save( CACHE_DIR + cached_file_name, stacked_data )
        self.load_file_timestamp_indeces.shape = ( len(self.load_file_timestamp_indeces),   )
        self.load_file_timestamps.shape = ( len(self.load_file_timestamps),   )


    # Retunrs 'True' or 'False'
    def load_file(self,datetime_stru):

        new_file_name = self.file_name_get(datetime_stru)

        if new_file_name == self.load_file_name:
            # Required file is currently loaded
            result = (hasattr(self,"load_file_data")) and (len(self.load_file_data) > 0)
            if (not result) and (not new_file_name in self.reported_missing_files):
                self.reported_missing_files += [new_file_name ]
                print( "{}file {} was not found, or was empty".format(self.info_lines_prefix, new_file_name))
            return result

        if new_file_name in self.not_found_file_names:
            # Required file was not found in a previous attempt
            return False

        if (self.load_file_name != "") and hasattr(self,"load_file_data") and (len(self.load_file_data) > 0):
            del self.load_file_data
            if self.load_file_timestamp_indeces is not None: del self.load_file_timestamp_indeces
            del self.load_file_timestamps
            self.load_file_data = np.array(()) # We create an empty array for empty / missing files  'len(self.load_file_data)) will return zero if no record isload'

        self.load_file_name = new_file_name
        self.load_file_name_base = os.path.basename(self.load_file_name)
        self.cursor_rec_index = 0

        if self.detailed_file_open_line:
            if self.data_descr == "":
                data_descr = "data from file "
            else:
                data_descr = f'"{self.data_descr}" data from file ' 
        else:
            data_descr = ""
        print("{}Reading {}\"{}\"...".format( self.info_lines_prefix, data_descr,  new_file_name), end='')

        fname, ext = os.path.splitext(new_file_name)
        if self.file_format == "":
            # Auto
            if (ext.lower()==".ssc") or (ext.lower()==".scc") :
                self.file_format = "SC"
            elif ext.lower()==".cdf" :
                self.file_format = "CDF"
            else:
                self.file_format = "ASCII-dat"

        assert self.file_format in SUPPORTED_FILE_FORMATS

        # Reading in file records and determining 'timeindeces' and timelabels
        if not os.path.isdir( self.path ):
            errmsg = "{}PATH \"{}\" DOES NOT EXIST!".format(self.info_lines_prefix, self.path)
            print(errmsg)
            self.not_found_file_names += [ new_file_name ]
            assert False, errmsg
            return False

        if not os.path.isfile( self.path + new_file_name):
            print(self.info_lines_prefix  + "FILE NOT FOUND!")
            self.not_found_file_names += [ new_file_name ]
            self.reported_missing_files += [ new_file_name ]
            return False

        store_in_buffer = True
        file_is_read_now = True

        # Reusing already read data, if we have it:
        if new_file_name in self.buffered_files.keys():
            self.load_file_data = self.buffered_files[new_file_name]["data"]
            self.load_file_timestamp_indeces = self.buffered_files[new_file_name]["timeindeces"]
            self.load_file_timestamps = self.buffered_files[new_file_name]["timestamps"]
            store_in_buffer = False
            file_is_read_now = False
            if len(self.load_file_data) > 0: 
                print("re-using buffered data with {} records... ".format(len(self.load_file_data)) , end='')                
            else:
                print("file was empty")

        elif (CACHE_DIR is not None):
            cached_file_name = self.path.replace("\\", "_" ).replace("/", "_" ).replace(':','_') + new_file_name+ "." +CACHE_FILE_EXT
            if os.path.isfile( CACHE_DIR + cached_file_name ) and (os.path.getmtime(CACHE_DIR + cached_file_name) > os.path.getmtime(self.path + new_file_name) ) :
                cache_stacked_data = np.load(CACHE_DIR + cached_file_name)
                if cache_stacked_data.shape[1] == self.nof_comps + 2:
                    print(" using cache file ({}) ".format(cached_file_name),end='')
                    # stacked_data = np.load(CACHE_DIR + cached_file_name)
                    self.load_file_timestamp_indeces = cache_stacked_data[...,0]  # Note: timeindeces will be rounded to integers later
                    self.load_file_timestamps  = cache_stacked_data[...,1]
                    self.load_file_data = cache_stacked_data[...,2:]
                    assert self.load_file_data.shape[1] == self.nof_comps,"Cached file has different number of components!"
                    file_is_read_now = False

        if file_is_read_now:
            # First-time reading of this file (or need other components that are in cache file)
            if (self.file_format == "SC"):
                # Reading compressed SC file
                errmsg, warnings, self.load_file_timestamps, self.load_file_data, comps = sc_file_reader.load(self.path + new_file_name, False, self.data_comps_file_colindeces)
                assert errmsg == "", f"Error reading SC file: {errmsg}"
                assert False, "Revise this part"
                for warning in warnings:
                    mu.warning(f"Reading SC file: {warning}")
                if self.file_record_interval_type == "msec":
                    self.load_file_timestamp_indeces = (self.load_file_timestamps % 1.0) * 86400. * 1000.
                elif self.file_record_interval_type == "sec":
                    self.load_file_timestamp_indeces = (self.load_file_timestamps % 1.0) * 86400.
                elif self.file_record_interval_type == "min":
                    self.load_file_timestamp_indeces = (self.load_file_timestamps % 1.0) * 1440.
                else:
                    assert False,"timeindeces calculation not yet implemented"

            elif (self.file_format == "IAGA-2002"):
                assert '#' in self.file_record_comments
                # dtype = [ "|S10", "|S12" ]  # Date and time fields
                dtype = [ "|S10", datetime.datetime ]  # Date-start and fully parsed time (combined with file date)
                for i in range(self.nof_comps): dtype += [ float ]   # Magnetic components
                
                time_format_mask = ["%H:%M:%S.000", "%H:%M:%S.%f"] [0] # Two options for IAGA files. For one second and one minute IAGA-2002 files, probably using the first is better
                time_converter = { 1: lambda s:datetime.datetime.combine(datetime_stru, datetime.datetime.strptime(s, time_format_mask).time()) }

                # For using numyp.genfromtxt() to load IAGA-2002 file, we need to know the number of optional comment line
                nr_of_obligatrory_header_lines = 13
                nr_of_user_comments = 0
                with open(self.path + new_file_name, encoding="utf-8") as f:
                    for line_nr, line in enumerate(f):
                        if line_nr < nr_of_obligatrory_header_lines: continue
                        if line[0:2] == " #":
                            nr_of_user_comments += 1
                        else:
                            # We should be in the compoment header
                            assert line[0:4] == "DATE"
                            break
                IAGA_2002_data = np.genfromtxt( self.path + new_file_name, 
                                                 skip_header=nr_of_obligatrory_header_lines + nr_of_user_comments, 
                                                 dtype = dtype, usecols=[0,1] + self.data_comps_file_colindeces, 
                                                 names = True, 
                                                 comments=self.file_record_comments, 
                                                 converters = time_converter, 
                                                 encoding='utf-8', 
                                                 missing_values = ["99999.00", "88888.00"],
                                                 ndmin = 1)
                if (len(IAGA_2002_data)>0) and (IAGA_2002_data[0][0]==b'DATE'):
                    # Data header (after file header and comment records) is read as a record:
                    IAGA_2002_data = IAGA_2002_data[1:]
                datetime_from_filename_str = datetime_stru.strftime("%Y-%m-%d")
                datetime_from_filename_bytes = bytes(datetime_from_filename_str, 'utf-8' )
                
                # Checking if all date string corresponds to file names date:
                assert np.all(  IAGA_2002_data[ IAGA_2002_data.dtype.names[0] ]  == datetime_from_filename_bytes), "Date mismatch in IAGA-2002 file"
                    
                # n = len(IAGA_2002_data)
                # self.load_file_data = np.empty( ( n, self.nof_comps )  )
                # self.load_file_timestamps = np.empty( n )
                self.load_file_timestamps = mu.datetime_2_exceltime( IAGA_2002_data[ IAGA_2002_data.dtype.names[1] ]) 
                
                # Transforming the components into a 2D array:
                mag_components = list(IAGA_2002_data.dtype.names[2:])  # skip first two fields
                assert len(mag_components) == len( self.data_comps_file_colindeces )
                assert len(mag_components) in [ 3, 4] # Only reading of three vector component, or reading of three vector component plus F is supported
                ordered_mag_components = [ "?" ] * len(mag_components)
                for mag_comp in mag_components:
                    # Format of component name: <IAGA-CODE><COMPNAME>
                    assert len(mag_comp) == 4
                    # Ordering data column in XYZ or HDZ order
                    if mag_comp[3] in "XH": 
                        ordered_mag_components[0] = mag_comp
                    elif mag_comp[3] in "EDY":  
                        ordered_mag_components[1] = mag_comp
                    elif mag_comp[3] in "Z":  
                        ordered_mag_components[2] = mag_comp
                    elif mag_comp[3] in "F":  
                        assert len(mag_components)  == 4
                        ordered_mag_components[3] = mag_comp
                    else:
                        assert False, 'Unknown magnetic component "{}" in IAGA-2002 file "{}"'.format(mag_comp, new_file_name)
               
                IAGA_2002_data_packed = rfn.repack_fields(IAGA_2002_data[ordered_mag_components])
                self.load_file_data = rfn.structured_to_unstructured(IAGA_2002_data_packed)
                    
                if self.file_record_interval_sec == 1:
                    self.load_file_timestamp_indeces = np.round((self.load_file_timestamps % 1.0) * 86400.).astype(int)
                elif self.file_record_interval_sec == 60:
                    self.load_file_timestamp_indeces = np.round((self.load_file_timestamps % 1.0) * 1440.).astype(int)
                else:
                    assert False,"timeindeces calculation not yet implemented"
                # first_comp_col_index = 0

            elif (self.file_format == "CDF"):
                import cdflib
                #  CDF file for NAQ
                cdf_file = cdflib.CDF(self.path + new_file_name)
                if self.file_record_interval_type == "sec":
                    cdf_time = cdf_file.varget('time')
                    N = 86400
                    assert len(cdf_time)==N
                    assert self.nof_comps == 5
                    H=cdf_file.varget('HNvar')
                    D=cdf_file.varget('HEvar')
                    Z=cdf_file.varget('Zvar')
                    T1=cdf_file.varget('T1')
                    T2=cdf_file.varget('T2')
                    HNscv = cdf_file.varget('HNscv')
                    HEscv = cdf_file.varget('HEscv')
                    Zscv  = cdf_file.varget('Zscv')
                    assert (len(HNscv) == len(HEscv) == len(Zscv))
                    for i in range(len(HNscv)):
                        assert HNscv[i]==HNscv[0]
                        assert HEscv[i]==HEscv[0]
                        assert Zscv[i]==Zscv[0]
                    assert (len(H)==N) and (len(D)==N) and (len(Z)==N) and (len(T1)==N)  and (len(T2)==N)
                    self.load_file_data = np.empty( (N, self.nof_comps)  )
                    self.load_file_data[:,0] = H*HNscv[0]
                    self.load_file_data[:,1] = D*HEscv[0]
                    self.load_file_data[:,2] = Z*Zscv[0]
                    self.load_file_data[:,3] = T1
                    self.load_file_data[:,4] = T2
                    datetime_from_filename = mu.datetime_2_exceltime( self.truncate_datetime_to_filename_period( datetime_stru )  )
                    self.load_file_timestamps = np.linspace(datetime_from_filename, datetime_from_filename+1.0, N, endpoint = False)
                    self.load_file_timestamp_indeces = (self.load_file_timestamps % 1.0) * 86400.

                    self.load_file_datas = np.empty( (N, self.nof_comps+1)  )
                    self.load_file_datas[:,0] = self.load_file_timestamps
                    self.load_file_datas[:,1] = H*HNscv[0]
                    self.load_file_datas[:,2] = D*HEscv[0]
                    self.load_file_datas[:,3] = Z*Zscv[0]
                    self.load_file_datas[:,4] = T1
                    self.load_file_datas[:,5] = T2
                    # np.savetxt("d:\\{}.csv".format(new_file_name), self.load_file_datas, delimiter=";", fmt=['%.8f', '%.3f', '%.3f', '%.3f', '%.3f' ,'%.3f'], header="datetime;H;D;Z;T1;T2")
                else:
                    assert False


            else:
                
               
                # Notes for 'genfromtxt()'
                # If 'dtype' parameter is used and it is a list/tuple, and/or 'names = True' is used, then 'np.genfromtxt()' function
                # returns and 1-dimensional array, where each element is a list (a record in file).
                # But we prefer using 2D array, where each row is a file record, each column is a component.
                # Using "ndim=2" parameter, maybe a 2D array is returned, but 'ndim' is supported only from 1.23.0 of numpy for 'genfromtxt()'.
                # Original implementation: 'dtype =  [int] * first_comp_col_index  + [float] * self.nof_comps' s
                
                # ToDo: read %A timestamp component as float
                
                genfromtxt_dtype = []
                for time_comp_index in range(len(self.time_components_file_colindeces)):
                    genfromtxt_dtype += [  (self.time_components_names[time_comp_index], 'i4') ]
                for data_comp_index in range(len(self.data_comps_file_colindeces)):
                    genfromtxt_dtype += [  (self.comp_names[data_comp_index], 'f8' ) ]
                
                # genfromtxt_dtype = ['i4'] * len(self.time_components_file_colindeces) + ['f8'] * len(self.data_comps_file_colindeces)
                
                self.load_file_data = np.genfromtxt( self.path + new_file_name, skip_header=self.fix_header_lines, dtype = genfromtxt_dtype, usecols=self.time_components_file_colindeces + self.data_comps_file_colindeces, comments=self.file_record_comments, ndmin = 1)
                datetime_from_filename = mu.datetime_2_exceltime( self.truncate_datetime_to_filename_period( datetime_stru )  )

                if self.file_record_interval_unit == "milliseconds":
                    assert self.time_components_file_colindeces[0] == self.FILE_COLINDEX_HOUR
                    assert self.time_components_file_colindeces[1] == self.FILE_COLINDEX_MIN
                    assert self.time_components_file_colindeces[2] == self.FILE_COLINDEX_SEC
                    assert self.time_components_file_colindeces[3] == self.FILE_COLINDEX_MSEC
                    self.load_file_timestamps        = self.calc_timestamps(datetime_from_filename, self.load_file_data[TIME_COMPONENT_NAME__HOUR], self.load_file_data[TIME_COMPONENT_NAME__MINUTE], self.load_file_data[TIME_COMPONENT_NAME__SECOND], self.load_file_data[TIME_COMPONENT_NAME__MILLISECOND])
                    self.load_file_timestamp_indeces = self.calc_timestamp_indeces(self.load_file_data[TIME_COMPONENT_NAME__HOUR], self.load_file_data[TIME_COMPONENT_NAME__MINUTE], self.load_file_data[TIME_COMPONENT_NAME__SECOND], self.load_file_data[TIME_COMPONENT_NAME__MILLISECOND])
                elif self.file_record_interval_unit == "seconds":
                    assert self.time_components_file_colindeces[0] == self.FILE_COLINDEX_HOUR
                    assert self.time_components_file_colindeces[1] == self.FILE_COLINDEX_MIN
                    assert self.time_components_file_colindeces[2] == self.FILE_COLINDEX_SEC
                    assert self.FILE_COLINDEX_MSEC is None
                    self.load_file_timestamps        = self.calc_timestamps(datetime_from_filename, self.load_file_data[TIME_COMPONENT_NAME__HOUR], self.load_file_data[TIME_COMPONENT_NAME__MINUTE], self.load_file_data[TIME_COMPONENT_NAME__SECOND])
                    self.load_file_timestamp_indeces = self.calc_timestamp_indeces(self.load_file_data[TIME_COMPONENT_NAME__HOUR], self.load_file_data[TIME_COMPONENT_NAME__MINUTE], self.load_file_data[TIME_COMPONENT_NAME__SECOND])
                elif self.file_record_interval_unit == "minutes":
                    assert self.time_components_file_colindeces[0] == self.FILE_COLINDEX_HOUR
                    assert self.time_components_file_colindeces[1] == self.FILE_COLINDEX_MIN
                    assert self.FILE_COLINDEX_SEC is None
                    assert self.FILE_COLINDEX_MSEC is None
                    self.load_file_timestamps        = self.calc_timestamps(datetime_from_filename, self.load_file_data[TIME_COMPONENT_NAME__HOUR], self.load_file_data[TIME_COMPONENT_NAME__MINUTE], 0)
                    self.load_file_timestamp_indeces = self.calc_timestamp_indeces(self.load_file_data[TIME_COMPONENT_NAME__HOUR], self.load_file_data[TIME_COMPONENT_NAME__MINUTE], 0)
                else:
                    assert False,"time-indeces calculation not yet implemented"
                    

                # Removing not used time label components ("hh mm ss") from 'self.load_file_data'
                # self.load_file_data = np.delete(self.load_file_data, range(first_comp_col_index), axis=1)
                # del first_comp_col_index
                
                self.load_file_data = rfn.drop_fields( self.load_file_data, self.time_components_names, usemask = False  )
                
                # ToDo: keep using structured array
                self.load_file_data =  rfn.structured_to_unstructured(self.load_file_data)

            if len(self.load_file_data)==0:
                print("file has no any record!")
                # Alwawy adding empty file to buffer:
                assert store_in_buffer
                self.buffered_files[new_file_name] = { "size": 0, "data": self.load_file_data, "timeindeces":self.load_file_timestamp_indeces, "timestamps": self.load_file_timestamps }
                return False

            # Saving not yet records-filtered file data (self.load_file_timestamps / self.load_file_timestamp_indeces / self.load_file_data) to cache files
            if CACHE_DIR is not None:
                self.save_cached_file(cached_file_name)

        self.nr_of_read_files += 1

        # Applying 'nan' and range selection filters for raw data
        
        #  Possibilites:
        #    1. In a record, one component value is invalid (e.g. 99999.99)
        #    1. In an ohter record, an other component value is invalid (e.g. 99999.99)
        #    2. In a record, two or more, but not all, records are invalud
        #    3. In a record all records are invalid

        # In the list range definition format, there should be a separate limit for all components.
        # Converting the list into a dict of components:
        for min_max in ['min','max']:
            if isinstance(self.range_filter[min_max], (int, float)):
                # One simple number for all commponents
                self.range_filter[min_max] = dict.fromkeys(self.comp_names, self.range_filter[min_max])
            elif type(self.range_filter[min_max])==list:
                self.range_filter[min_max] = dict( zip(self.comp_names, self.range_filter[min_max])  )
        
        # We search the special case, if all outliers is a given value, typically 99999.99
        range_out_single_value__per_comp = dict.fromkeys( self.comp_names, None)  # Will be set to "99999.99", of the np.nan, if there are more range out components
        range_out_single_value__all_comps = None  # Will be np.nan, if there are more different range out components
        
        nr_invalid_componts__per_comp   =  dict.fromkeys( self.comp_names, 0)
        nr_invalid_componts__per_record =  np.zeros( len(self.load_file_data) )
        
        nr_of_read_records = len(self.load_file_data)
        
        for comp_index, comp_name in enumerate(self.comp_names):
            comp_data = self.load_file_data[:, comp_index]
            
            # Checking nans:
            nan_record_indeces = np.isnan(comp_data)
            # Checking min/max:
            out_of_range_record_indeces = np.full( nr_of_read_records, False )
            for min_max in ['min','max']:
                range_limit = self.range_filter[min_max].get(comp_name, None)
                if range_limit:
                    if min_max == "min":
                        out_of_range_record_indeces |= (comp_data < range_limit)
                    elif min_max == "max":
                        out_of_range_record_indeces |= (comp_data > range_limit)
                    else:
                        assert False
                    # Note: 'nans' will not be included in 'out_of_range_record_indeces'
            if np.any(out_of_range_record_indeces):
                all_invalid_record_indeces = nan_record_indeces | out_of_range_record_indeces
                range_out_values = comp_data[out_of_range_record_indeces]
                if np.all( range_out_values == range_out_values[0] ):
                    # E.g. "99999.99"
                    range_out_single_value__per_comp[comp_name] = range_out_values[0]
                    if range_out_single_value__all_comps is None:
                        # First set of 'range_out_single_value__all_comps'
                        range_out_single_value__all_comps = range_out_single_value__per_comp[comp_name]
                    elif (not np.isnan(range_out_single_value__all_comps)) and (range_out_single_value__all_comps != range_out_single_value__per_comp[comp_name]):
                        # Different range out values at different compoments
                        range_out_single_value__all_comps = np.nan
                else:
                    range_out_single_value__per_comp[comp_name] = np.nan
                    range_out_single_value__all_comps = np.nan
                if self.leave_partial_records:
                    # Inavlidating out of range values, as the record may contain other components with valid values, so records will be kept
                    comp_data[out_of_range_record_indeces] = np.nan
            else:
                all_invalid_record_indeces = nan_record_indeces 
            if np.any(all_invalid_record_indeces):
                nr_invalid_componts__per_comp[comp_name] = np.count_nonzero(all_invalid_record_indeces)
                nr_invalid_componts__per_record[all_invalid_record_indeces] += 1
            
        all_comps_are_invalid_record_indeces = (nr_invalid_componts__per_record == self.nof_comps)
        
        nr_records_to_delete  = np.count_nonzero(all_comps_are_invalid_record_indeces)
        assert self.leave_partial_records
        nr_partial_records    = np.count_nonzero( nr_invalid_componts__per_record[ ~all_comps_are_invalid_record_indeces ] > 0 )
        
        if nr_records_to_delete > 0:
            assert file_is_read_now  # Data from RAM buffer may already be filtered!
            if nr_records_to_delete == nr_of_read_records:
                msg = f'all the {nr_of_read_records} read records are empty/invalid => not using the file'
            else:
                msg = f'file-records:{nr_of_read_records} empty-records:{nr_records_to_delete}'
            self.load_file_data               = np.delete(self.load_file_data, all_comps_are_invalid_record_indeces, axis=0)                
            self.load_file_timestamps         = np.delete(self.load_file_timestamps, all_comps_are_invalid_record_indeces)
            if self.load_file_timestamp_indeces is not None: # It is None in case of not aligned time stamps
                self.load_file_timestamp_indeces  = np.delete(self.load_file_timestamp_indeces, all_comps_are_invalid_record_indeces)

            assert (len(self.load_file_data) == 0) == ( nr_records_to_delete == nr_of_read_records)
        else:
            msg = f'file-records:{nr_of_read_records}'
        if (nr_records_to_delete < nr_of_read_records) and (nr_partial_records>0):
            msg += f' partial-records:{nr_partial_records}'
            
        
        
        if (range_out_single_value__all_comps is not None) and not(np.isnan(range_out_single_value__all_comps)):
            msg += f' (nan:"{range_out_single_value__all_comps}")'
            
        there_is_valid_record = len(self.load_file_timestamps)  >  0
        if not there_is_valid_record:
            # All records are deleted!            
            msg+= " - no valid record remained!"


        if file_is_read_now: print(msg, end='')
        
        if not there_is_valid_record:
            assert store_in_buffer
            self.buffered_files[new_file_name] = { "size": 0, "data": self.load_file_data, "timeindeces":self.load_file_timestamp_indeces, "timestamps": self.load_file_timestamps }
            return False
            
        if store_in_buffer:  # and (len(self.load_file_data)>0 )
            # Storing even the empty files (to not to read the empty records again)
            # A new file will be read -> buffering currently load file data:
            assert not new_file_name in self.buffered_files.keys()
            required_new_buff_size = self.load_file_data.nbytes + self.load_file_timestamps.nbytes + ( 0 if self.load_file_timestamp_indeces is None else self.load_file_timestamp_indeces.nbytes ) 
            while self.sum_buffer_data_size + required_new_buff_size > MAX_BUFFERED_SIZE:
                oldest_buffered_file, oldest_buffered_content = self.buffered_files.popitem(last=False) 
                print("------- FILE RAM BUFFER FULL, DELETING {} FILE DATA  ---------".format(oldest_buffered_file))
                self.sum_buffer_data_size -= oldest_buffered_content["size"]
                # self.buffered_files.pop(  self.buffered_files[filedata_to_delete] ) 

            self.buffered_files[new_file_name] = \
                {
                  "size"       : required_new_buff_size,
                  "data"       : self.load_file_data,
                  "timeindeces": self.load_file_timestamp_indeces,
                  "timestamps" : self.load_file_timestamps
                }
                
            self.sum_buffer_data_size +=   required_new_buff_size
            

        # if len(self.AR_raw_min)>0:
        #     assert len(self.AR_raw_min) == self.nof_comps
        #     # orig_nof_records = len(self.load_file_data)
        #     delete_record_indeces = np.full( len(self.load_file_data), False )  # Creates a boolean array filled with "False" values
        #     for icomp in range(self.nof_comps):
        #          deleted_record_indeces_of_comp = (self.load_file_data[:, icomp] < self.AR_raw_min[icomp])
        #          if np.any(deleted_record_indeces_of_comp):
        #              delete_record_indeces = delete_record_indeces  | deleted_record_indeces_of_comp
        #     if delete_record_indeces.any(): # There is at least one record to delete
        #         self.load_file_timestamps   = np.delete(self.load_file_timestamps, delete_record_indeces)
        #         if self.load_file_timestamp_indeces is not None:
        #             self.load_file_timestamp_indeces  = np.delete(self.load_file_timestamp_indeces, delete_record_indeces)
        #         self.load_file_data         = np.delete(self.load_file_data, delete_record_indeces, axis=0)
        #         nof_removed_records += np.count_nonzero(delete_record_indeces)
        #         applied_filters = applied_filters + ["AR_raw_min"]

        # if len(self.AR_raw_max)>0:
        #     assert len(self.AR_raw_max) == self.nof_comps, "len(AR_raw_max)={:} nof_comps={:}".format(len(self.AR_raw_max), self.nof_comps)
        #     # orig_nof_records = len(self.load_file_data)
        #     delete_record_indeces = np.full( len(self.load_file_data), False )  # Creates a boolean array filled with "False" values
        #     for icomp in range(self.nof_comps):
        #         delete_record_indeces = delete_record_indeces  | (self.load_file_data[:, icomp] > self.AR_raw_max[icomp])
        #     if delete_record_indeces.any(): # There is at least one record to delete
        #         self.load_file_timestamps   = np.delete(self.load_file_timestamps, delete_record_indeces)
        #         if self.load_file_timestamp_indeces is not None:
        #             self.load_file_timestamp_indeces  = np.delete(self.load_file_timestamp_indeces, delete_record_indeces)
        #         self.load_file_data         = np.delete(self.load_file_data, delete_record_indeces, axis=0)
        #         nof_removed_records += np.count_nonzero(delete_record_indeces)
        #         applied_filters = applied_filters + ["AR_raw_max"]

        #if nof_removed_records > 0:
        #    print("  {:} record{} removed due to filter(s) ""{}""".format( nof_removed_records, mu.plural_suffix(nof_removed_records, "s"), applied_filters), end='')
        #    if len(self.load_file_data)==0:
        #        print(" no valid records remained!")
        #        return False

    	# Timeindeces must be integer values, so rounding it to int:
        if self.load_file_timestamp_indeces is not None:
            assert isinstance( self.load_file_timestamp_indeces[0], (int, np.uint32, np.int64) ) 

        assert min(self.load_file_timestamps) > 44000
        assert max(self.load_file_timestamps) < 56000
        
        if self.load_file_timestamp_indeces is not None:
            assert len(self.load_file_timestamp_indeces) == len(self.load_file_timestamps)
            # Checking if timelabels (time indeces) are monolithically increasing:
            # First, checking for "same-time" time indeces
            d_timeindeces = np.diff( self.load_file_timestamp_indeces )
            zero_d_timeindeces = (d_timeindeces == 0)  # Return a True/False array
            if warn_zero_timeindeces_diff and np.any(zero_d_timeindeces):
                print(self.info_lines_prefix + "There is/are similar time-labels/time-indeces:")
                for i, value in enumerate(zero_d_timeindeces):
                    if value:
                        print("timelabel:{:.5f}  timeindex:{:d}".format(self.load_file_timestamps[i], round(self.load_file_timestamp_indeces[i])) )
                mu.warning("")

            accepted_time_jumps = (d_timeindeces >= 0)
            if not np.all(accepted_time_jumps):
                print(self.info_lines_prefix+ "There is/are decreasing time labels!")
                print(self.info_lines_prefix+ "Records with time back jumps: ")
                not_accepted_time_jump_indeces = [i for i in range(len(accepted_time_jumps)) if not accepted_time_jumps[i]]
                for i in not_accepted_time_jump_indeces:
                    print(self.info_lines_prefix  + "timelabel:{:.5f}   index-jump:{:d} sec".format(self.load_file_timestamps[i], d_timeindeces[i]) )
                assert False
        else:
            pass
            # ToDo: check if 'self.load_file_timestamps' is monolitically increasing

        # Applying conversion rules:
        if self.conversion_defs != None:
            for conversion_def in self.conversion_defs:
                assert len(conversion_def) == 2
                if conversion_def[0]=="M":
                    assert isinstance(conversion_def[1],(float,int))
                    self.load_file_data =   self.load_file_data * conversion_def[1]
                elif conversion_def[0]=="O":  # Offset
                    assert isinstance(conversion_def[1],(float,int))
                    self.load_file_data =   self.load_file_data + conversion_def[1]
                else:
                    assert False, "Unknown conversion: {}".format(conversion_def[0])

        print()
        return True

    # Internally called when no data found for the use requested method.
    # It returns the same staff for all cases of 'no record found' situation.
    def no_data_found(self, for_time, case_id, in_case_of_missing_data, extra_info = None):
        if in_case_of_missing_data == "drop-error":
            if isinstance(for_time,(list,tuple)):
                for_time = for_time[0]
            if not isinstance(for_time,datetime.datetime):
                assert isinstance(for_time, float)
                for_time = mu.exceltime_to_datetime(for_time, "sec")
            if isinstance(for_time, datetime.datetime):
                dt = for_time
            else:
                dt = mu.exceltime_to_datetime(for_time)
            assert False, "No record found in ""{}"" for date/time {:04d}.{:02d}.{:02d} {:02d}:{:02d}:{:02d} (case_id:{})".format(self.load_file_name_base, dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, case_id )
        if in_case_of_missing_data == "drop-warning":
            assert False, "not yet implemented"
        if in_case_of_missing_data in [ "none"  ]: # ,  "just-return-nans"
            errmsg = "No valid data could be found" if self.nr_of_read_files>0 else "No data file found"         
            return errmsg, np.nan, self.no_record_found_return_data, [], extra_info
        assert False, f"{in_case_of_missing_data=}"

    def get_record_at_cursor(self):
        return self.load_file_data[ self.cursor_rec_index ]

    # Measning of 'raw_data_result_timestamp_diff': raw (file) data is reported, if file and request timelabel's difference is up to 'raw_data_result_timestamp_diff'    
    def set_processing_parameters(self, raw_data_result_timestamp_diff = "auto", interpolated_data_result_halfwindow_size = None,  min_samples_for_avg_calc = DEF_MIN_SAMPLES_FOR_AVG_CALC):
        self.min_samples_for_avg_calc = DEF_MIN_SAMPLES_FOR_AVG_CALC 
        # If 'interpolated_data_result_halfwindow_size' is None, then it is not set at this call!
        if interpolated_data_result_halfwindow_size is not None: self.interpolated_data_result_halfwindow_size = interpolated_data_result_halfwindow_size  # if 'None', then interpolation is not used!        
        if raw_data_result_timestamp_diff == "auto":
            if self.interpolated_data_result_halfwindow_size is None:
                # Interpolation is not used, so reporting raw data only
                self.raw_data_result_timestamp_diff = self.timeindex_interval_in_sec  / 86400.0 * 0.501
            else:
                self.raw_data_result_timestamp_diff =  self.timeindex_interval_in_sec  / 86400.0 * 0.1
        else:
            self.raw_data_result_timestamp_diff = raw_data_result_timestamp_diff
            

        
        

    # Main function of 'dat_file_reader_class': it reads out data from data files for the desired timestamps/time labels.
    # 'time_label_def' must be either a single float, an array/list/tuple of floats, or sometimes a range, defined as a list.
    # Returns four variables:
    # 1. ErrorMsg: None if there is at least one record with one valid value. If all requested records components are invalid, then a related message
    # 2. timelabels as float array 
    # 3. Data as 2-D np.ndarry with  'self.nof_comps' 2nd axis  (filled with NaN, if no data found)
    # 4. List (in a numpy array) of non complete record indeces
    # 5. Optionally: extra information: generally it is 'None'. Used for 'calc_stat'.
    # Possible values for 'in_case_of_missing_data': 'drop-error', 'just-return-nans'
    # 'master_call' is True when this function is called from outside, then it is False in second and third recursice call. It is used for logging of missing files.
    def get_data(self, time_label_def, result_data_type = file_read_result_data_type.records_for_timestamps,  in_case_of_missing_data = "drop-error" , calc_stat=False, master_call = True, dtype = None):

        assert isinstance(result_data_type, file_read_result_data_type) 

        assert in_case_of_missing_data in ["drop-error", "drop-warning", "none" ]

        if master_call:
            self.reported_missing_files = []
            # Checking if directory exist:

            current_path = self.path        
            if isinstance(time_label_def, float):
                first_timelabel = time_label_def
            else:
                first_timelabel = time_label_def[0]
            
            if self.path_has_datetime_mask: current_path = mu.exceltime_2_datetime(first_timelabel).strftime(current_path)
            assert os.path.isdir(current_path), "Data path \"{}\" not found".format(current_path)
            assert not self.path_has_datetime_mask, "Having mask in the path part is not yet supported"
            self.nr_of_read_files = 0

        # if isinstance(time_label_def, (list,tuple,np.ndarray)):
        #     if ( (result_data_type != file_read_result_data_type.arithm_avg_within_period) and (result_data_type != file_read_result_data_type.all_within_period) ) or isinstance(time_label_def[0], (list,tuple)):
        #         # More timelabels (or intervals) are supplied, we read and collect data for each of them
        #         for i, i_time_label_def in enumerate(time_label_def):
        #             timelabel, data, extra_data = self.get_data(i_time_label_def, file_read_result_data_type.records_for_timestamps, in_case_of_missing_data, calc_stat, master_call = False)
        #             result_timelabels[i] = timelabel
        #             result_data[i][:] = data[:]
        #             # ToDo: collect extra-data
        #             if dtype is not None:
        #                 # Converting 2D array into structured array:
        #                 # HIBÁS / TODO: ez a blokk törött (nem 'dtype is not None' hívással sosem tesztelt).
        #                 # 'result_data' egy sima 2D tömb, 'result_data[comp]' string-indexelés hibát dob rá;
        #                 # 'result_data_as_structured_array[:,i]' szintén hibás indexelés egy 1D struktúrált tömbön;
        #                 # és az irány is fordítva van (a friss, üres struktúrált tömbből másolna, nem oda).
        #                 # JAVÍTANDÓ, mielőtt bárki 'dtype' paraméterrel hívja a get_data()-t!
        #                 assert len(dtype)==self.nof_comps
        #                 result_data_as_structured_array = np.empty( len(timelabel), dtype = dtype )
        #                 for i,comp in enumerate(dtype):
        #                     result_data[comp][:] = result_data_as_structured_array[:,i]
        #                 result_data = result_data_as_structured_array

        #         return result_timelabels, result_data, None if not calc_stat else extra_data
            
        match result_data_type:
            
          case file_read_result_data_type.records_for_timestamps:
            # --------------------------------------------------------------------
            # Most common call type: caller requests record(s) for given time stamp(s)
            # Trying to used record exactly for 'time_label'
            # 'time_label_def' can be a simple timestamp, of a list/tuple/array of timestamps 
            # --------------------------------------------------------------------
            if not isinstance(time_label_def, (list,tuple,np.ndarray)):
                assert isinstance(time_label_def, float)
                time_label_def = [ time_label_def ]
            n = len(time_label_def)

            result_data       = np.full( (n, self.nof_comps), np.nan )
            # result_timelabels = np.empty( n )            
            # result_timelabels[:] = time_label_def[:]
            result_timelabels = time_label_def.copy()
            there_is_valid_result = False
            
            
            if self.interpolated_data_result_halfwindow_size:
                round_to_interval_sec = None # Using interpolation, so not rounding to fixed times
            else:
                round_to_interval_sec = self.timeindex_interval_in_sec
            
            datetime_timestamps = mu.exceltime_2_datetime(time_label_def, round_to_interval_sec)
            
            for result_data_index, time_label in enumerate(time_label_def):
                
                # timelabel, data, extra_data = self.get_data(i_time_label_def, file_read_result_data_type.records_for_timestamps, in_case_of_missing_data, calc_stat, master_call = False)
                datetime_timestamp = datetime_timestamps[result_data_index]
                if not self.load_file( datetime_timestamp ): 
                    continue
                
                # First step: trying to use original (not interpolated) data
                using_aligned_timestamps = False
                if (self.file_record_interval_aligned) and (datetime_timestamp.microsecond == 0) and (type(self.timeindex_interval_in_sec) is int):
                    
                    if (datetime_timestamp.second == 0):
                        # ToDo: rvise this par
                        using_aligned_timestamps = self.timeindex_interval_in_sec < 60
                    else:
                        using_aligned_timestamps =  (datetime_timestamp.second % self.timeindex_interval_in_sec) == 0
                    
                if using_aligned_timestamps:                     
                    
                    # Searching record based on time stamp index 
                    requested_timeindex = self.calc_timestamp_indeces(datetime_timestamp.hour, datetime_timestamp.minute, datetime_timestamp.second)
                    if requested_timeindex in self.load_file_timestamp_indeces:
                        load_data_index = self.load_file_timestamp_indeces.searchsorted(requested_timeindex)
                        assert abs (result_timelabels[result_data_index] - self.load_file_timestamps[load_data_index]) * 86400 < 0.1
                        result_data[result_data_index,:] = self.load_file_data[load_data_index,:]
                else:
                    closest_record_index = find_closest_index(self.load_file_timestamps, time_label)
                    dtime = abs(self.load_file_timestamps[closest_record_index] - time_label)
                    if dtime <= self.raw_data_result_timestamp_diff:
                        result_data[result_data_index,:] = self.load_file_data[closest_record_index,:]
                        
                if not np.any( np.isnan(result_data[result_data_index]) ):
                    # All components could be read
                    there_is_valid_result = True
                else:
                    # Trying to fill missing compoments by interpolation
                    if (not self.interpolated_data_result_halfwindow_size is None):
                        self.set_record_with_time_interpolation(time_label, result_data[result_data_index])
                    if not np.all( np.isnan(result_data[result_data_index]) ):
                        there_is_valid_result = True # There is at least one valid components
                    
            rows_with_at_least_one_nan = np.isnan(result_data).any(axis=1)  # use "numpy.all" to get records with full nan's 
            noncomplete_record_indeces = np.where(rows_with_at_least_one_nan)[0]
                        
            if dtype is not None:
                # Converting 2D array into structured array
                # Not tested!
                assert len(dtype)==self.nof_comps
                result_data_as_structured_array = np.empty( len(result_timelabels), dtype = dtype )
                for comp_index,comp_def in enumerate(dtype):
                    result_data_as_structured_array[comp_def[0]][:] = result_data[:,comp_index] 
                result_data = result_data_as_structured_array

            if not there_is_valid_result:
                assert len(noncomplete_record_indeces)  == len(time_label_def) 
                errmsg = "No valid data could be found" if self.nr_of_read_files > 0 else "No data file found"         
            else:
                errmsg = None
            if not calc_stat: extra_data = None
            return errmsg, result_timelabels, result_data, noncomplete_record_indeces, extra_data
        
            
          case file_read_result_data_type._1min_arithm_avg:
            assert isinstance(time_label_def,float)
            timelabels, data_in_avg_window = self.get_all_data_in_time_window( time_label_def - 0.5 / 1440.0,  time_label_def + 0.5 / 1440.0, min_samples = self.min_samples_for_avg_calc)
            if len(timelabels)==0:
                return self.no_data_found(time_label_def, "no-data-in-windows", in_case_of_missing_data)
            return np.mean(timelabels), np.mean(data_in_avg_window, axis = 0), None

          case  file_read_result_data_type.arithm_avg_within_period | file_read_result_data_type.all_within_period:
            assert len(time_label_def) == 2
            if isinstance(time_label_def[0], datetime.datetime): time_label_def[0]=mu.datetime_2_exceltime(time_label_def[0])
            if isinstance(time_label_def[1], datetime.datetime): time_label_def[1]=mu.datetime_2_exceltime(time_label_def[1])
            assert time_label_def[1] > time_label_def[0]
            timelabels, data_in_avg_window = self.get_all_data_in_time_window(time_label_def[0],  time_label_def[1], min_samples = self.min_samples_for_avg_calc )
            if len(timelabels)==0:
                return self.no_data_found(time_label_def, "no-data-in-windows", in_case_of_missing_data,  [0, np.nan, np.nan] if calc_stat else None)
            if calc_stat:
                # Calculating and reporting standard deviation and slope:
                stddev = np.std(data_in_avg_window, axis = 0)
                slope =  np.polyfit(timelabels,data_in_avg_window,1)[0]
                # Slope is in <data_unit>/day, converting it to <data_unit>/sec:
                slope =  slope / 86400.0
                # res = scipy.stats.linregress(timelables, data_in_avg_window)
                extra_data = [len(timelabels),  stddev, slope]
            else:
                extra_data = None
            if (result_data_type == file_read_result_data_type.all_within_period):
                return timelabels, data_in_avg_window, extra_data
            else:
                return np.mean(timelabels), np.mean(data_in_avg_window, axis = 0), extra_data


        assert False


    # Fills emppty data for the given timelabel by interpolation of the two neighbouring records.
    # Typically called from 'get_data()'.
    # For interpolation, we use a 'window' that is +-self.interpolated_data_result_halfwindow_size around the time label.
    # 'record' must be an 1-D array, it is an input and output parameter as well.
    # Only 'nan' values in 'record' will be set.
    def set_record_with_time_interpolation(self, time_label, record):
        
        comp_indeces_to_set = np.isnan(record)
        comp_indeces_to_set = np.where( comp_indeces_to_set  )[0]
        
        timelabels_in_window, data_in_window = self.get_all_data_in_time_window( time_label - self.interpolated_data_result_halfwindow_size , time_label + self.interpolated_data_result_halfwindow_size, min_samples = 1)
        n_window = len(data_in_window)

        if n_window==0: return
            # return self.no_data_found(time_label,"any-data-in-window",in_case_of_missing_data)

        if n_window == 1:
            # Can not apply interpolation
            record[comp_indeces_to_set] = data_in_window[0][comp_indeces_to_set]
            return
        
        
        closest_index = np.searchsorted(timelabels_in_window, time_label, side =  "left") # find_closest_index(timelabels_in_window, time_label
        if closest_index == n_window: closest_index = n_window-1

        for comp_index in comp_indeces_to_set:
            
            if timelabels_in_window[-1] <= time_label:
                assert closest_index == n_window-1
                # All of the data found in the interpolation window is at or before the timelabel. We can not interpolat, just use a value
                for j in range( n_window-1, -1, -1):
                    if not np.isnan( data_in_window[j,comp_index] ):
                        record[comp_index] = data_in_window[j,comp_index]
                        break
            elif timelabels_in_window[0] >= time_label:
                # All of the data found in the window is after the timelabel. 
                assert closest_index == 0
                for j in range( 0, n_window, +1):
                    if not np.isnan( data_in_window[j,comp_index] ):
                        record[comp_index] = data_in_window[j,comp_index]
                        break
            else: 
                assert timelabels_in_window[closest_index-1] <= time_label
                assert timelabels_in_window[closest_index  ] >= time_label
                # Values before/after can be nan's
                before_index_non_nan, after_index_non_nan = None, None
                for before_index in range( max(0,closest_index-1),-1,-1):
                    if not np.isnan(data_in_window[before_index][comp_index]): 
                        before_index_non_nan = before_index
                        break
                for after_index in range( closest_index,n_window,+1):
                    #print(after_index,comp_index)
                    #if after_index == 9:
                    #    pass
                    if not np.isnan(data_in_window[after_index][comp_index]): 
                        after_index_non_nan = after_index
                        break
                if before_index_non_nan is None:
                    if after_index_non_nan is not None:
                        record[comp_index] = data_in_window[after_index_non_nan, comp_index]
                elif after_index_non_nan is None:
                    record[comp_index] = data_in_window[before_index_non_nan, comp_index]
                else:
                    record[comp_index] = np.interp( time_label, [timelabels_in_window[before_index_non_nan], timelabels_in_window[after_index_non_nan]],  [ data_in_window[before_index_non_nan,comp_index],  data_in_window[after_index_non_nan,comp_index] ] )

        return # We return the results in 'record'


    # Returns all data within a window.
    # Returns: timelabels as 1D 'horizontal' array, and data as 2D array, where each row is a record, each column is a component
    # If not data found, then it returns empty arrays. It does not drop error, if there is no data!
    def get_all_data_in_time_window(self, datetime_from, datetime_to, min_samples = 1):

        if isinstance(datetime_from, datetime.datetime): datetime_from = mu.datetime_2_exceltime(datetime_from)
        if isinstance(datetime_to, datetime.datetime):   datetime_to   = mu.datetime_2_exceltime(datetime_to)

        # Increasing the window by half of the record interval
        # Note: this should be done only if datetime_from / datetime_to resolution is equal to file record interval!
        # Todo: consider comparing using "timeindeces" method - but it is tricky if requtested window passes file period

        datetime_from = datetime_from - self.timeindex_interval_in_sec / 2. / 86400.
        datetime_to   = datetime_to   + self.timeindex_interval_in_sec / 2. / 86400.

        # Step-1: searching the first record at the start of the time window
        dt_from       = mu.exceltime_2_datetime(datetime_from, self.timeindex_interval_in_sec)
        dt_to         = mu.exceltime_2_datetime(datetime_to,   self.timeindex_interval_in_sec)

        # Determining dates for file names
        # Note: 'mu.day_range()' could also be used?!
        if self.file_name_get(dt_from) == self.file_name_get(dt_to):
            read_file_dates = [ dt_from ]
        else:
            read_file_dates = [ dt_from ]
            assert self.file_name_interval_type == "days"
            dt_next_day = dt_from + datetime.timedelta(days=1)
            while True:
              read_file_dates += [ dt_next_day ]
              if dt_next_day>=dt_to: break
              dt_next_day = dt_next_day + datetime.timedelta(days=1)

        joined_data = None
        joined_timelabels = None

        for read_file_date in read_file_dates:

            if not self.load_file(read_file_date):
                continue

            first_record_index_at_or_after_the_start_time = np.searchsorted(self.load_file_timestamps,datetime_from, "left")
            last_record_index_at_or_before_end_time       = np.searchsorted(self.load_file_timestamps,datetime_to,   "right")

            assert first_record_index_at_or_after_the_start_time <= last_record_index_at_or_before_end_time

            if first_record_index_at_or_after_the_start_time==len(self.load_file_data):
                # All read records are before the time window low side
                continue

            if last_record_index_at_or_before_end_time == 0:
                # 'datetime_to' is before the first record of the file: all records are after the high side of the window
                continue


            if first_record_index_at_or_after_the_start_time>0:
                assert( self.load_file_timestamps[first_record_index_at_or_after_the_start_time-1] < datetime_from)
            assert( self.load_file_timestamps[first_record_index_at_or_after_the_start_time ] >= datetime_from)

            # Checking of the applied algorithms:
            assert(  self.load_file_timestamps[first_record_index_at_or_after_the_start_time]  >= datetime_from )
            if first_record_index_at_or_after_the_start_time>0:
                assert(  self.load_file_timestamps[first_record_index_at_or_after_the_start_time-1]  < datetime_from )


            if first_record_index_at_or_after_the_start_time == last_record_index_at_or_before_end_time:
                # No data within window situation
                assert not (datetime_from <= self.load_file_timestamps[first_record_index_at_or_after_the_start_time] <=datetime_to)
                continue

            last_record_index_at_or_before_end_time = last_record_index_at_or_before_end_time - 1

            # if last_record_index_at_or_before_end_time == len(self.load_file_data):
            #     last_record_index_at_or_before_end_time  - 1

            assert(  self.load_file_timestamps[last_record_index_at_or_before_end_time]  <= datetime_to )
            if last_record_index_at_or_before_end_time < len(self.load_file_data)-1:
                assert(  self.load_file_timestamps[last_record_index_at_or_before_end_time + 1]  > datetime_to )

            if joined_data is None:
                joined_timelabels = self.load_file_timestamps[ first_record_index_at_or_after_the_start_time:last_record_index_at_or_before_end_time + 1 ]
                joined_data       = self.load_file_data[ first_record_index_at_or_after_the_start_time:last_record_index_at_or_before_end_time  + 1]
            else:
                joined_timelabels = np.append(joined_timelabels, self.load_file_timestamps[ first_record_index_at_or_after_the_start_time:last_record_index_at_or_before_end_time+ 1 ])
                joined_data       = np.append(joined_data,       self.load_file_data[ first_record_index_at_or_after_the_start_time:last_record_index_at_or_before_end_time +1 ], axis = 0  )

            if len(joined_timelabels) != len(joined_data):
                print (f"{first_record_index_at_or_after_the_start_time=}  {last_record_index_at_or_before_end_time=}")
                print (f"{len(self.load_file_data)=}")
                print (f"{len(self.load_file_timestamps)=}")
                print(self.load_file_data[ first_record_index_at_or_after_the_start_time:last_record_index_at_or_before_end_time  + 1])
                assert len(joined_timelabels) == len(joined_data), f"{len(joined_timelabels)=}/{len(joined_data)=}"
            if  len(joined_timelabels) == 0:
                assert  False
            assert np.all( joined_timelabels[:] >= datetime_from)
            assert np.all( joined_timelabels[:] <= datetime_to)

        if joined_data is None:
            return np.empty([0]), np.empty([0, self.nof_comps])
        else:
            return joined_timelabels, joined_data




# Returns one-minute interval time labels between two dates.
def get_minute_timelabels(time_start, time_end):
    assert isinstance(time_start,float)
    # Rounding 'time_start' and 'time_end' to nearest minutes

    # Nem pontos elso probalkozas:

    # time_start = round(time_start * 1440) / 1440
    # time_end   = round(time_end   * 1440) / 1440
    # nof_minutes = round ((time_end - time_start)*1440 )
    # return np.linspace(time_start, time_end, num = nof_minutes)

    dt_start = mu.exceltime_2_datetime(time_start,60)
    assert dt_start.second == 0
    assert dt_start.microsecond == 0

    dt_end = mu.exceltime_2_datetime(time_end, 60)
    assert dt_end.second == 0
    assert dt_end.microsecond == 0

    # if dt_end.second > 30:
    #    dt_end = dt_end + dt.timedelta(minutes=1)
    # dt_end.replace(second=0, microsecond = 0)

    n0 = int  ((time_end - time_start) * 1440.0) + 2
    result = np.empty( n0 )

    n=0

    while (dt_start <= dt_end) :
        if n>=n0:
            assert False
        result[n] = mu.datetime_2_exceltime(dt_start)
        dt_start = dt_start + datetime.timedelta(minutes=1)
        # Just to be sure that dt_start is rounded to minutes:
        dt_start.replace(second=0, microsecond = 0)
        n += 1

    result = result[:n]
    return result


def get_second_timelabels(time_start, time_end):
    assert isinstance(time_start,float)
    # Rounding 'time_start' and 'time_end' to nearest minutes

    # Nem pontos elso probalkozas:

    # time_start = round(time_start * 1440) / 1440
    # time_end   = round(time_end   * 1440) / 1440
    # nof_minutes = round ((time_end - time_start)*1440 )
    # return np.linspace(time_start, time_end, num = nof_minutes)

    dt_start = mu.exceltime_to_datetime(time_start, "sec")
    assert dt_start.second == 0
    assert dt_start.microsecond == 0

    dt_end = mu.exceltime_to_datetime(time_end, "sec")
    assert dt_end.second == 0
    assert dt_end.microsecond == 0

    # if dt_end.second > 30:
    #    dt_end = dt_end + dt.timedelta(minutes=1)
    # dt_end.replace(second=0, microsecond = 0)

    n0 = int  ((time_end - time_start) * 86400.0) + 2
    result = np.empty( n0 )

    n=0

    while (dt_start <= dt_end) :
        if n>=n0:
            assert False
        result[n] = mu.datetime_2_exceltime(dt_start)
        dt_start = dt_start + datetime.timedelta(seconds=1)
        # Just to be sure that dt_start is rounded to minutes:
        # dt_start.replace(second=0, microsecond = 0)
        n += 1

    result = result[:n]
    return result


# Creates an instance of 'dat_file_reader_class' from a config definition, 
# that is typically read from YAML config file.
# Supported formats:
#  - string, with comma separated fields: "<patf>,<file-name-mask>,... (see source)
#  - dictionary   
def dat_file_reader_class_from_config_def(config_def):

    if type(config_def) is str:
        fields =  config_def.split()
        path = fields[0]
        file_name_mask = fields[1]
        file_interval_def = fields[2]
        components = mu.ensure_list(fields[3])
        colindeces = mu.ensure_list(fields[4])
        assert False, "Not fully implemented"        
    elif type(config_def) is dict:
        path = config_def["path"]
        file_name_mask = config_def["file-name-mask"]
        file_interval_def = config_def["interval"]
        components = config_def.get("components")
        colindeces = config_def.get("column-indeces")
        hour_colindex = config_def.get("hour-colindex", None)
        range_filter_def = config_def.get("range", None)
        
    else:
        assert False, "{} format based definicion is not supported".format(type(config_def))
    res = dat_file_reader_class(path, file_name_mask, file_interval_def, components, colindeces, FILE_COLINDEX_HOUR = hour_colindex, range_filter = range_filter_def  )
    return res

