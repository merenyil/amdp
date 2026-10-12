# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 08:45:05 2021

@author: Merényi László
"""

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.dates as mdates
import numpy as np
import misc_utils as mu
import matplotlib.cm as cm  # for 'rainbow' color
import matplotlib
import datetime

normal_figsize = (5, 4)  # width / height
small_figsize = (5, 2)   # width / height

ticks_for_deg_axis =  np.arange(0,360+90,90)

def rainbow_color_get(index, nof_all_items):
    if nof_all_items == 1:
        return 'black'
    else:
        return cm.rainbow( index / (nof_all_items-1) )



# 'Y' is an 2D-array/list/tuple, containing one or more components!
# 'label' is a list/tuple: label for each components - it can be empty
# 'fmt' can be a list or a simple string varaiable
# 'vertical_separators' can be '%years'
def plot_datetime_y(timestamps, Y, Yaxis_label, fmt='-', datetimeformatter = "%Y-%m-%d %H:%M:%S", color='blue', label = None, b_close_plot = True, b_configure_datetime_axis = True, remove_mean = False, markersize = None, vertical_separators = None):

    # python_datetime = np.copy(_datetime)
    series_names = label
    if not isinstance(timestamps, np.ndarray): timestamps = np.array(timestamps)

    # if not  np.all( _datetime[~np.isnan(_datetime)] ) > 44560.0:
    #     for i in range(len(_datetime)):
    #         if np.isnan(_datetime[i]):
    #             print(i,"nan")
    #         elif _datetime[i] < 44560.0:
    #             print(_datetime[i])
    #     assert np.all( _datetime > 44560.0)
    # assert np.all( _datetime[~np.isnan(_datetime)] < 46000.0)

    python_datetime = mu.exceltime_2_matplotlib_time(timestamps)

    Y = mu.ensure_list(Y)

    for i,y in enumerate(Y):
        assert(len(y) == len(timestamps))
        _fmt= fmt[i] if isinstance( fmt, (list, tuple)) else fmt
        _color = color[i] if isinstance( color, (list, tuple)) else color
        _series_names = series_names[i] if isinstance( series_names, (list, tuple)) else (series_names if series_names is not None else "")
        if remove_mean:
            plt.plot(python_datetime, y - np.mean(y), _fmt, color=_color, label = _series_names, markersize = markersize)
        else:
            plt.plot(python_datetime, y, _fmt, color=_color, label = series_names, markersize = markersize)
    plt.xlabel("date/time")
    plt.ylabel(Yaxis_label)
    plt.ticklabel_format(axis='y', style='plain', scilimits=None, useOffset=False)

    if b_configure_datetime_axis:

        ax = plt.gca()   #     get axes
        [ xmin, xmax ] =  ax.get_xlim()

        # date_range_day = max(python_datetime) - min(python_datetime)
        date_range_day = xmax - xmin

        formatter = mdates.DateFormatter(datetimeformatter)
        ax.xaxis.set_major_formatter(formatter)

        nr_of_xaxis_ticks = 3.0
        if date_range_day < 0.05:
            interval = round ( (date_range_day * 86400) / nr_of_xaxis_ticks )  # 'nr_of_xaxis_ticks' pieces of locators for "yyyy:mm:dd hh:mm" format
            locator = mdates.SecondLocator( interval = interval ) #     set hour as tick marks
            ax.xaxis.set_major_locator(locator)
        elif date_range_day < 0.5:
            interval = round ( (date_range_day * 1440) / nr_of_xaxis_ticks )  # 'nr_of_xaxis_ticks' pieces of locators for "yyyy:mm:dd hh:mm" format
            locator = mdates.MinuteLocator( interval = interval ) #     set hour as tick marks
            ax.xaxis.set_major_locator(locator)
        elif date_range_day < 4:
            interval = round ( (date_range_day * 24) / nr_of_xaxis_ticks)  # 'nr_of_xaxis_ticks' pieces of locators for "yyyy:mm:dd hh:mm" format
            locator = mdates.HourLocator( interval = interval ) #     set hour as tick marks
            ax.xaxis.set_major_locator(locator)
        else:
            interval = round(  date_range_day / nr_of_xaxis_ticks )
            locator = mdates.DayLocator(interval = interval)
            ax.xaxis.set_major_locator(locator)


        if isinstance(label, (list, tuple)) and (len(label)>1): plt.legend()
        
    if (vertical_separators is not None) and (type(vertical_separators) is str):
        if vertical_separators == "%years":
            min_year = mu.exceltime_2_datetime(min(timestamps)).year
            max_year = mu.exceltime_2_datetime(max(timestamps)).year
            running_year = datetime.datetime(min_year,1,1)
            max_year = datetime.datetime(max_year,1,1)
            vertical_separators = []
            while running_year<=max_year:
                vertical_separators += [ mu.datetime_2_exceltime(running_year) ]
                running_year = running_year.replace(year = running_year.year + 1)
        
    if (vertical_separators is not None) and (len(vertical_separators)>0):
        ax = plt.gca()
        vertical_separators = mu.exceltime_2_matplotlib_time(vertical_separators)
            
        ax.xaxis.set_minor_locator( matplotlib.ticker.FixedLocator(vertical_separators) )
        plt.grid(visible=True, which="minor", axis='x', linestyle='dotted' )
        

    if b_close_plot:

        plt.show()


# 'X' is a list or 1D array. Or just an integer value: in that case it is turned into range(X).
# 'Y' can be a  simple l1D array, or a list of 1D array
# symbol/color can also be a list
# Other way to define groups is to use 'point_groups'
def plot_xy(X, Y, Xaxis_label, Yaxis_label, symbol='.', color='blue', Xaxis_range = None, Yaxis_range = None, xaxis_datetimeformatter = "", close_plot = True, point_groups = None, label = None, title = None, linewidth = 1, vertical_separators = None, plot_width = None, horizontal_separators = None):
    if isinstance(X,int): X=range(X)
    
    series_names = label

    if plot_width is not None:
        # fig = plt.figure(figsize=(FIGURE_SIZE[0], FIGURE_SIZE[1]*3), layout = ['constrained', 'compressed', 'tight', 'none', None][2] )
        fig = plt.figure(figsize=(plot_width, 4.8) )
        
    b_more_series = isinstance(Y, (list, tuple)) and isinstance(Y[0], (list, tuple, np.ndarray)) 

    if point_groups == None:
        if b_more_series:
            # Multi-component plot:
            for i,y in enumerate(Y):
                if isinstance(symbol,list):
                    act_symbol = symbol[i]
                else:
                    act_symbol = symbol
                if isinstance(color, list):
                    act_color = color[i]
                else:
                    act_color = color
                if series_names is not None:
                    act_label = series_names[i]
                else:
                    act_label = Yaxis_label
                plt.plot(X, y, act_symbol, color = act_color, label=act_label, linewidth = linewidth)
                if title is not None:
                    plt.title(title)
        else:
            # Simple one component plot
            plt.plot(X, Y, symbol, color = color, linewidth = linewidth);
        plt.xlabel(Xaxis_label)
        plt.ylabel(Yaxis_label)
        plt.ticklabel_format(axis='both', style='plain', scilimits=None, useOffset=False)
        if isinstance(series_names, (list, tuple)): plt.legend(framealpha=1, frameon=True, bbox_to_anchor=(1.1, 1.05))
    else:
        assert len(point_groups) == 3  # See below
        group_indeces = point_groups[0]  # 'group_indeces' should be a list (or range), where each items are the indeces, that belongs to one group
        group_labels  = point_groups[1]
        group_symbols = point_groups[2]
        assert len(group_indeces) == len(group_labels)
        assert len(group_indeces) == len(group_symbols)
        for i_group in range( len(group_indeces)  ):
            plt.plot( [ X[i] for i in group_indeces[i_group] ], [ Y[j] for j in group_indeces[i_group]], group_symbols[i_group], label = group_labels[i_group] )
        plt.xlabel(Xaxis_label)
        plt.ylabel(Yaxis_label)
        plt.legend(framealpha=1, frameon=True, bbox_to_anchor=(1.1, 1.05))

    if Xaxis_range != None:
        axes = plt.gca()
        axes.set_xlim([Xaxis_range[0],Xaxis_range[1]])
    if Yaxis_range != None:
        axes = plt.gca()
        axes.set_ylim([Yaxis_range[0],Yaxis_range[1]])
    if xaxis_datetimeformatter != "":
        assert False, "This part is still required? 'plot_datetime_y()' could be used instead?!"  # Added on 2025.03.23
        ax = plt.gca()   #     get axes
        formatter = mdates.DateFormatter(xaxis_datetimeformatter) #     format as dates
        ax.xaxis.set_major_formatter(formatter)
        if xaxis_datetimeformatter == "%H:%M":
            locator = mdates.HourLocator() #     set hour as tick marks
            ax.xaxis.set_major_locator(locator)
        else :
            locator = mdates.DayLocator()
            ax.xaxis.set_major_locator(locator)

    if title is not None:
        plt.title(title)

    if (vertical_separators is not None) and (len(vertical_separators)>0):
        ax = plt.gca()
        ax.xaxis.set_minor_locator( matplotlib.ticker.FixedLocator(vertical_separators) )
        plt.grid(visible=True, which="minor", axis='x', linestyle='dotted' )

    if (horizontal_separators is not None) and (len(horizontal_separators)>0):
        # Implemented in 2025.05.10, but not working
        ax = plt.gca()
        ax.yaxis.set_minor_locator( matplotlib.ticker.FixedLocator(horizontal_separators) )
        plt.grid(visible=True, which="minor", axis='y', linestyle='dotted' )

    if close_plot:
        plt.show()


def plot_size_set(size):
    if (size=='small'):
        plt.figure(figsize=small_figsize)
    elif (size=='normal'):
        plt.figure(figsize=normal_figsize)
    else:
        assert False

# Source data originally lie on a evently spaced grid, so we do not need to use interpolation.
# Used for grid type absolute measurements.
# Measurements must not be in X/Y increasing order
# But, after ordering X/Y: 'X and Y must both be ordered monotonically'.
def contour_plot_from_gridded_data(title, x,y,z, x_div_count, y_div_count, ilevels = 100):
    assert x_div_count * y_div_count == len(x)
    assert len(x) == len(y)
    assert len(x) == len(z)

    x_points = np.empty(x_div_count)
    y_points = np.empty(y_div_count)
    z_points = np.empty((x_div_count,y_div_count))
    z_points[:] = np.nan

    # First, collecting out x grid points - in unordered manner
    x_points[0] = min(x)
    collected_x_count = 1
    for i in range(len(x)):
        this_gridpoint_exist = False
        for k  in range(collected_x_count):
            if x_points[k]==x[i]:
                this_gridpoint_exist = True
                break
        if not this_gridpoint_exist:
            assert(collected_x_count < x_div_count)
            x_points[collected_x_count] = x[i]
            collected_x_count = collected_x_count + 1
    assert( collected_x_count == x_div_count)

    # Then, collecting out y grid points - in unordered manner
    y_points[0] = min(y)
    collected_y_count = 1
    for i in range(len(y)):
        this_gridpoint_exist = False
        for k  in range(collected_y_count):
            if y_points[k]==y[i]:
                this_gridpoint_exist = True
                break
        if not this_gridpoint_exist:
            assert(collected_y_count < y_div_count)
            y_points[collected_y_count] = y[i]
            collected_y_count = collected_y_count + 1
    assert( collected_y_count == y_div_count)

    # Then sorting x and y grid points
    x_points = np.sort(x_points)
    y_points = np.sort(y_points)

    # Then collecting z points in the grid order:
    for i in range(len(z)):
        for ix in range(x_div_count):
            if x[i] == x_points[ix]: break
        for iy in range(y_div_count):
            if y[i] == y_points[iy]: break
        assert(np.isnan(z_points[iy,ix]))
        z_points[iy,ix]=z[i]

    # grid_X, grid_Y = np.meshgrid(x_points, y_points)
    color_map = plt.get_cmap('seismic')
    plt.contourf(x_points, y_points, z_points, levels = ilevels, cmap =  color_map)
    plt.colorbar()
    # plt.xlim([0,360])
    # plt.ylim([0,360])
    # plt.pause(0.5)
    plt.title(title)
    plt.show()



def contour_plot(x,y,z, x_label, y_label, z_title, x_minmax, y_minmax, ngridx=100, ngridy=100, xticks = None, yticks = None):

    # Methods:
    # 1) 'https://matplotlib.org/2.0.2/examples/pylab_examples/contour_demo.html'
    # 2) Using Kriging: https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/examples/01_universal.html#sphx-glr-examples-01-universal-py
    # 3) 'https://matplotlib.org/stable/gallery/images_contours_and_fields/irregulardatagrid.html'

    # This draws a 3D surface with one color:
    # hf = plt.figure()
    # ha = hf.add_subplot(111, projection='3d')
    # ha.plot_trisurf(x, y, z)
    # plt.xlabel('hc (°)')
    # plt.ylabel('vc (°)')
    # plt.show()


    n = len(x)
    x_min = x_minmax[0]
    x_max = x_minmax[1]
    y_min = y_minmax[0]
    y_max = y_minmax[1]

    LEVELS = 20
    b_plot_with_triang = False
    b_plot_with_tricontour = True

    # Grid values:
    xi = np.linspace(x_min, x_max, ngridx)
    yi = np.linspace(y_min, y_max, ngridy)

    if b_plot_with_triang:
        # 3/1) Using 'tri.Triangulation()'
        # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
        fig, (ax1) = plt.subplots(nrows=1)
        triang = tri.Triangulation(x, y)
        interpolator = tri.LinearTriInterpolator(triang, z)
        Xi, Yi = np.meshgrid(xi, yi)
        zi = interpolator(Xi, Yi)  # 'zi' will have the 'z' values interpolated to the grid

        ax1.contour(xi, yi, zi, levels=LEVELS, linewidths=0.5, colors='k')
        cntr1 = ax1.contourf(xi, yi, zi, levels=14, cmap="RdBu_r")

        fig.colorbar(cntr1, ax=ax1)
        ax1.plot(x, y, 'k.', ms=1)
        ax1.set(xlim=(x_min, x_max), ylim=(y_min, y_max))
        ax1.set_title("{} with triang".format(z_title))


    if b_plot_with_tricontour:
        # 3/2) - directly using 'tricontour()'/'tricontourf()' built in functions.
        fig = plt.figure( figsize=(10,10) )
        ax2 = plt.subplot(1,1,1)
        ax2.tricontour(x, y, z, levels=LEVELS, linewidths=0.5, colors='k')
        cntr2 = ax2.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

        fig.colorbar(cntr2, ax=ax2)
        ax2.plot(x, y, 'ko', ms=3)
        ax2.set(xlim=(x_min, x_max), ylim=(y_min, y_max))
        ax2.set_title('{} with tricontour'.format(z_title))

        # plt.subplots_adjust(hspace=0.5)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if not xticks is None: plt.xticks(xticks)
    if not yticks is None: plt.yticks(yticks)
    plt.grid()
    plt.show()



