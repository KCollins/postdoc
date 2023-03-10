# functions for visualization of magnetometer data. 

# Importing packages
import pandas as pd
smag = __import__('supermag-api')          # SuperMAG python API
logon = 'kd8oxt'                              # SuperMAG ID

import plotly.express as px                # for mapping, mainly
import plotly.graph_objects as go          # for mapping, mainly

import matplotlib.pyplot as plt
import numpy as np

# For pulling data from CDAweb:
from ai import cdas
import datetime
from matplotlib import pyplot as plt

# For saving images:
import os

############################################################################################################################### 
# # MAGFIG Function to create stacked plot of conjugate magnetometers.
#
#
# INPUTS:
#
#    parameter = The parameter of interest - Bx, By, or Bz. North/South, East/West, and vertical, respectively.
#    start, end = datetimes of the start and end of plots
#
# OUTPUTS:
#
#    Figure of stacked plots for date in question, with events marked.


def magfig(
    parameter = 'Bz',
    start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
    end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
    maglist_a = ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],
    maglist_b = ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],
    is_displayed = False,
    is_saved = False, 
    events=None, event_fontdict = {'size':20,'weight':'bold'}
):
        # Magnetometer parameter dict so that we don't have to type the full string: 
        d = {'Bx':'MAGNETIC_NORTH_-_H', 'By':'MAGNETIC_EAST_-_E','Bz':'VERTICAL_DOWN_-_Z'}
        if is_saved:
            fname = 'output/' +str(start) + '_' +  str(parameter) + '.png'
            if os.path.exists(fname):
                print('Looks like ' + fname + ' has already been generated.')
                return 
                # raise Exception('This file has already been generated.')
        fig, axs = plt.subplots(len(maglist_a), figsize=(25, 25), constrained_layout=True)
        print('Plotting data for ' + str(len(maglist_a)) + ' magnetometers: ' + str(start))
        for idx, magname in enumerate(maglist_a):   # Plot Arctic mags:
            print('Plotting data for Arctic magnetometer #' + str(idx+1) + ': ' + magname.upper())
            try:                
                data = cdas.get_data(
                    'sp_phys',
                    'THG_L2_MAG_'+ magname.upper(),
                    start,
                    end,
                    ['thg_mag_'+ magname]
                )
                data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
                x =data['UT']
                y =data[d[parameter]]
                y = reject_outliers(y) # Remove power cycling artifacts on, e.g., PG2.
                axs[idx].plot(x,y)#x, y)
                axs[idx].set(xlabel='Time', ylabel=magname.upper())
                
                if events is not None:
                    # print('Plotting events...')
                    trans       = mpl.transforms.blended_transform_factory(axs[idx].transData,axs[idx].transAxes)
                    for event in events:
                        evt_dtime   = event.get('datetime')
                        evt_label   = event.get('label')
                        evt_color   = event.get('color','0.4')

                        axs[idx].axvline(evt_dtime,lw=1,ls='--',color=evt_color)
                        if evt_label is not None:
                            axs[idx].text(evt_dtime,0.01,evt_label,transform=trans,
                                    rotation=90,fontdict=event_fontdict,color=evt_color,
                                    va='bottom',ha='right')
                
                
                try: 
                    magname = maglist_b[idx]
                    ax2 = axs[idx].twinx()
                    print('Plotting data for Antarctic magnetometer #' + str(idx+1) + ': ' + magname.upper())
                    data = cdas.get_data(
                        'sp_phys',
                        'THG_L2_MAG_'+ magname.upper(),
                        start,
                        end,
                        ['thg_mag_'+ magname]
                    )
                    data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
                    x =data['UT']
                    y =data[d[parameter]]

                    color = 'tab:red'
                    # ax2.set_ylabel('Y2-axis', color = color)
                    # ax2.plot(y, dataset_2, color = color)
                    # ax2.tick_params(axis ='y', labelcolor = color)
                    y = reject_outliers(y) # Remove power cycling artifacts on, e.g., PG2.
                    ax2.plot(x,-y, color=color)#x, y)
                    ax2.set_ylabel(magname.upper(), color = color)
                    ax2.tick_params(axis ='y', labelcolor = color)
                except Exception as e:
                    print(e)
                    continue
            except Exception as e:
                print(e)
                continue
        fig.suptitle(str(start) + ' ' +  str(parameter), fontsize=30)    # Title the plot...
        if is_saved:
            print("Saving figure. " + fname)
            # fname = 'output/' +str(start) + '_' +  str(parameter) + '.png'
            fig.savefig(fname, dpi='figure', pad_inches=0.3)
        if is_displayed:
            return fig # TODO: Figure out how to suppress output here

###############################################################################################################################        
# Function to reject outliers. We'll need this to eliminate power cycling artifacts in the magnetometer plots.
def reject_outliers(y):   # y is the data in a 1D numpy array
    mean = np.mean(y)
    sd = np.std(y)
    final_list = np.copy(y)
    for n in range(len(y)):
        final_list[n] = y[n] if y[n] > mean - 3 * sd else np.nan
        final_list[n] = final_list[n] if final_list[n] < mean + 5 * sd else np.nan
    return final_list