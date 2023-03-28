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
import matplotlib as mpl

# For power plots:
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from scipy import signal
from scipy.fft import fft
from scipy.signal import butter, filtfilt, stft, spectrogram
import datetime as dt

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

############################################################################################################################### 
# # MAGSPECT Function to create power plots for conjugate magnetometers.
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

def magspect(
    parameter = 'Bx',
    start = datetime.datetime(2016,1,25,0,0), 
    end = datetime.datetime(2016,1,25,9,0), 
    maglist_a = ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],
    maglist_b = ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],
    is_displayed = False,
    is_saved = True, 
    events=None, event_fontdict = {'size':20,'weight':'bold'}, 
    myFmt = mdates.DateFormatter('%H:%M')
):
        # Magnetometer parameter dict so that we don't have to type the full string: 
        d = {'Bx':'MAGNETIC_NORTH_-_H', 'By':'MAGNETIC_EAST_-_E','Bz':'VERTICAL_DOWN_-_Z'}
        if is_saved:
            fname = 'output/' +str(start) + '_' +  str(parameter) + '.png'
            if os.path.exists(fname):
                print('Looks like ' + fname + ' has already been generated.')
                return 
                # raise Exception('This file has already been generated.')
        fig, axs = plt.subplots(len(maglist_a), 2, figsize=(25, 25), constrained_layout=True)
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
                y = reject_outliers(y)
                # Interpolate NaNs away: 
                df = pd.DataFrame(y, x)
                df = df.interpolate('linear')
                y = df[0].values
                # fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,10))
                xlim=[start,end]
                # ax1.plot(G13_dt,G13_Br)
                # ax1.plot(x, y-np.mean(y))
                # ax1.set_ylabel('[nT]') # CHECK: is this the correct unit?
                # ax1.set_xlabel('UT')
                # ax1.set(xlim=xlim)
                # myFmt = mdates.DateFormatter('%H:%M')
                axs[idx, 0].xaxis.set_major_formatter(myFmt)

                # f, t, Zxx = stft(G13_Br, fs=1/0.5, nperseg=1800,noverlap=1200)#
                # f, t, Zxx = stft(y-np.mean(y), fs=1, nperseg=1800,noverlap=1200)# Xueling's values
                f, t, Zxx = stft(y-np.mean(y), fs=1, nperseg = 256, noverlap = 128)#
                dt_list = [start+dt.timedelta(minutes=ii) for ii in t]

                cmap=axs[idx, 0].pcolormesh(dt_list, f*1000., np.abs(Zxx)*np.abs(Zxx), vmin=0, vmax=0.5)

                # ax2.axis([0, len(G13_Br)/2., 0, 40])
                #ax2.axis([0, len(y)/2., 0, 20])
                #plt.colorbar(cmap,label='spectrum [nT^2]')
                axs[idx, 0].set(xlim=xlim)
                axs[idx, 0].set(ylim=[0,20])
                axs[idx, 0].xaxis.set_major_formatter(myFmt)
                axs[idx, 0].set_title('STFT Power Spectrum: ' + magname.upper())
                axs[idx, 0].set_ylabel('Frequency [mHz]')
                axs[idx, 0].set_xlabel('UT')
                if events is not None:
                    # print('Plotting events...')
                    trans       = mpl.transforms.blended_transform_factory(axs[idx,0].transData,axs[idx,0].transAxes)
                    for event in events:
                        evt_dtime   = event.get('datetime')
                        evt_label   = event.get('label')
                        evt_color   = event.get('color','0.4')

                        axs[idx, 0].axvline(evt_dtime,lw=1,ls='--',color=evt_color)
                        if evt_label is not None:
                            axs[idx, 0].text(evt_dtime,0.01,evt_label,transform=trans,
                                    rotation=90,fontdict=event_fontdict,color=evt_color,
                                    va='bottom',ha='right')
            except Exception as e:
                print(e)
                continue
        for idx, magname in enumerate(maglist_b):   # Plot Antarctic mags:
            print('Plotting data for Antarctic magnetometer #' + str(idx+1) + ': ' + magname.upper())
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
                y = reject_outliers(y)
                # Interpolate NaNs away: 
                df = pd.DataFrame(y, x)
                df = df.interpolate('linear')
                y = df[0].values
                # fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,10))
                xlim=[start,end]
                # ax1.plot(G13_dt,G13_Br)
                # ax1.plot(x, y-np.mean(y))
                # ax1.set_ylabel('[nT]') # CHECK: is this the correct unit?
                # ax1.set_xlabel('UT')
                # ax1.set(xlim=xlim)
                # myFmt = mdates.DateFormatter('%H:%M')
                axs[idx, 1].xaxis.set_major_formatter(myFmt)

                # f, t, Zxx = stft(G13_Br, fs=1/0.5, nperseg=1800,noverlap=1200)#
                f, t, Zxx = stft(y-np.mean(y), fs=1, nperseg=1800,noverlap=1200)#
                dt_list = [start+dt.timedelta(seconds=ii) for ii in t]

                cmap=axs[idx, 1].pcolormesh(dt_list, f*1000., np.abs(Zxx)*np.abs(Zxx), vmin=0, vmax=0.5)

                # ax2.axis([0, len(G13_Br)/2., 0, 40])
                #ax2.axis([0, len(y)/2., 0, 20])
                #plt.colorbar(cmap,label='spectrum [nT^2]')
                axs[idx, 1].set(xlim=xlim)
                axs[idx, 1].set(ylim=[0,20])
                axs[idx, 1].xaxis.set_major_formatter(myFmt)
                axs[idx, 1].set_title('STFT Power Spectrum: ' + magname.upper())
                axs[idx, 1].set_ylabel('Frequency [mHz]')
                axs[idx, 1].set_xlabel('UT')
                if events is not None:
                    # print('Plotting events...')
                    trans       = mpl.transforms.blended_transform_factory(axs[idx,1].transData,axs[idx,1].transAxes)
                    for event in events:
                        evt_dtime   = event.get('datetime')
                        evt_label   = event.get('label')
                        evt_color   = event.get('color','0.4')

                        axs[idx, 1].axvline(evt_dtime,lw=1,ls='--',color=evt_color)
                        if evt_label is not None:
                            axs[idx, 1].text(evt_dtime,0.01,evt_label,transform=trans,
                                    rotation=90,fontdict=event_fontdict,color=evt_color,
                                    va='bottom',ha='right')
            except Exception as e:
                print(e)
                continue
        fig.suptitle(str(start) + ' ' +  str(parameter), fontsize=30)    # Title the plot...
        if is_saved:
            fname = 'output/PowerSpectrum_' +str(start) + '_' +  str(parameter) + '.png'
            print("Saving figure. " + fname)
            fig.savefig(fname, dpi='figure', pad_inches=0.3)
        if is_displayed:
            return fig # TODO: Figure out how to suppress output here
        
        
############################################################################################################################### 
# # MAGDF Function to create multi-indexable dataframe of all mag parameters for a given period of time. 
#
#
# INPUTS:
#
#    start, end = datetimes of the start and end of plots
#
# OUTPUTS:
#
#    Dataframe of Bx, By, Bz for each magnetometer in list.

def magdf(
    start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
    end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
    maglist_a = ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],  # Arctic magnetometers
    maglist_b = ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],  # Antarctic magnetometers
    is_saved = False, 
):
        # Magnetometer parameter dict so that we don't have to type the full string: 
        d = {'Bx':'MAGNETIC_NORTH_-_H', 'By':'MAGNETIC_EAST_-_E','Bz':'VERTICAL_DOWN_-_Z'}
        d_i = dict((v, k) for k, v in d.items()) # inverted mapping for col renaming later
        if is_saved:
            fname = 'output/' +str(start) + '_' + '.csv'
            if os.path.exists(fname):
                print('Looks like ' + fname + ' has already been generated.')
                return 
        UT = pd.date_range(start, end, freq ='S')    # preallocate time range
        full_df = pd.DataFrame(UT, columns=['UT'])   # preallocate dataframe
        full_df['UT'] = full_df['UT'].astype('datetime64[s]') # enforce 1s precision
        full_df['Magnetometer'] = ""
        for mags in [maglist_a, maglist_b]:
            for idx, magname in enumerate(mags):   # For each magnetometer, pull data and merge into full_df:
                print('Pulling data for magnetometer: ' + magname.upper())
                try:                
                    data = cdas.get_data(
                        'sp_phys',
                        'THG_L2_MAG_'+ magname.upper(),
                        start,
                        end,
                        ['thg_mag_'+ magname]
                    )
                    data['UT'] = pd.to_datetime(data['UT'])# unit='s')
                    df = pd.DataFrame(data)
                    df.rename(columns=d_i, inplace=True)    # mnemonic column names
                    
                    df['Magnetometer'] = magname.upper()
                    full_df = pd.concat([full_df, df])
                    
                    # print(df)
                except Exception as e:
                    print(e)
                    continue
        full_df['UT'] = full_df['UT'].astype('datetime64[s]') # enforce 1s precision
        full_df.drop(columns = ['UT_1']) # discard superfluous column
        full_df = full_df[full_df['Magnetometer'] != ''] # drop empty rows
        if is_saved:
            print('Saving as a CSV.')
            full_df.to_csv(fname)
        # print(full_df)
        return full_df 