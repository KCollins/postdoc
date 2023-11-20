# functions for visualization of magnetometer data. 

# Importing packages
import pandas as pd
smag = __import__('supermag-api')        # SuperMAG python API
logon = 'kd8oxt'                              # SuperMAG ID

import plotly.express as px    # for mapping, mainly
import plotly.graph_objects as go        # for mapping, mainly

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
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, hann
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# For saving images:
import os

# For fill_nan:
from scipy import interpolate

# for Shane's code to pull from TGO:
import requests

############################################################################################################################### 


def magfetch(
    start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
    end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
    magname = 'atu', 
    is_verbose = False, 
    tgopw = 'ResUseNoCom'
):
    """
    MAGFETCH 
        Function to fetch data for a given magnetometer. Pulls from ai.cdas or DTU.

        Arguments:
            start, end   : datetimes of the start and end of sampled data range.
            magname      : IAGA ID for magnetometer being sampled. e.g.: 'upn'
            is_verbose   : Boolean for whether debugging text is printed.
            tgopw        : Password for Tromsø Geophysical Observatory

        Returns:
            df           : pandas dataframe with columns ['UT', 'MAGNETIC_NORTH_-_H', 'MAGNETIC_EAST_-_E', 'VERTICAL_DOWN_-_Z']
    """
    # Pull password for TGO from local .txt file:
    file = open("tgopw.txt", "r")
    tgopw = file.read()
    if(is_verbose): print('Found Tromsø Geophysical Observatory password: ' + tgopw)
    if(magname in ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']): 
        if(is_verbose): print('Collecting data for ' + magname + ' from TGO.')
        if(end<start): print('End is after start... check inputs.')
        data = magfetchtgo(start, end, magname, tgopw = tgopw)
    else:
        data = cdas.get_data(
            'sp_phys',
            'THG_L2_MAG_'+ magname.upper(),
            start,
            end,
            ['thg_mag_'+ magname]
        )
    if(is_verbose): print('Data for ' + magname + ' collected: ' + str(len(data['UT'])) + ' samples.')
    return data

############################################################################################################################### 

def magfetchtgo(start, end, magname, tgopw = 'ResUseNoCom'):
    """
    Pulls data from a RESTful API with a link based on the date.

    Args:
        start (datetime.datetime): The start date of the data to be fetched.
        end (datetime.datetime): The end date of the data to be fetched.
        magname (str): The name of the magnetometer station.
        tgopw (str): Password for Tromso Geophysical Observatory.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing the fetched data.
    """

    df = pd.DataFrame()

    # Loop over each day from start to end
    for day in range(start.day, end.day + 1):
        # Generate the URL for the current day
        url = f'https://flux.phys.uit.no/cgi-bin/mkascii.cgi?site={magname}4d&year={start.year}&month={start.month}&day={day}&res=10sec&pwd='+ tgopw + '&format=XYZhtml&comps=DHZ&getdata=+Get+Data'

        # Fetch the data for the current day
        foo = pd.read_csv(url, skiprows = 6, delim_whitespace=True, usecols=range(5), index_col=False)
        # Convert the 'DD/MM/YYYY HH:MM:SS' column to datetime format
        foo['DD/MM/YYYY HH:MM:SS'] = foo['DD/MM/YYYY'] + ' ' + foo['HH:MM:SS']
        foo['UT'] = pd.to_datetime(foo['DD/MM/YYYY HH:MM:SS'], format='%d/%m/%Y %H:%M:%S')

        # Rename the columns
        foo.rename(columns={'X': 'MAGNETIC_NORTH_-_H', 'Y': 'MAGNETIC_EAST_-_E', 'Z': 'VERTICAL_DOWN_-_Z'}, inplace=True)
        
        df = pd.concat([df, foo])

    # # Convert the dataframe to a dictionary
    data = {
        'UT': df['UT'].to_numpy(),
        'MAGNETIC_NORTH_-_H': df['MAGNETIC_NORTH_-_H'].to_numpy(),
        'MAGNETIC_EAST_-_E': df['MAGNETIC_EAST_-_E'].to_numpy(),
        'VERTICAL_DOWN_-_Z': df['VERTICAL_DOWN_-_Z'].to_numpy()
    }
    # print(type(df))
    # return df
    return data
############################################################################################################################### 
def _download_remote_file(datetime: dt.datetime, site='tro', verbose=False):
    """Downloads file for a given date/site/orientation to the configured storage path

    Args:
        datetime (dt.datetime): date/time to download
        site (str, optional): Which site code
        orientation (str, optional): sensor orientation ('XYZ', 'HEZ')
    """
    def _build_request(datetime, site='tro'):
        request_dict = {#'site':config.site_codes[site],
                        'year':datetime.strftime('%Y'),
                        'month':'{}'.format(datetime.month),
                        'day':'{}'.format(datetime.day),
                        'res':'10sec',
                        'pwd': 'ResUseNoCom',
                        'format':'XYZhtml',
                        'comps':'SXYZ',
                        'getdata':' Get Data '}
        return request_dict
    
    # initiate request on remote server
    input_url = 'http://flux.phys.uit.no/cgi-bin/mkascii.cgi?'

    request_dict = _build_request(datetime, site=site)
    init_resp = requests.get(input_url, params=request_dict, timeout=(30,300))
    
    # check if header and footer are as expected
    if init_resp.text.startswith('<pre>\n') and init_resp.text.endswith('</pre>\n'):
        savedata = init_resp.text.lstrip('<pre>\n').rstrip('</pre>\n') + '\n'
        # find a place to put the data
        savefilepath = 'output/' + datetime.strftime(datetime) + '_' + site + '_' + sample_rate +'.txt'
        if(verbose): print(savefilepath)
        try:
            with open(savefilepath, 'w') as savefile:
                savefile.write(savedata)
        except FileNotFoundError:
            # os.makedirs(datetime.strftime(config.fpath_format.format(base=config.datapath_local)))
            with open(savefilepath, 'w') as savefile:
                savefile.write(savedata)
        return savedata
    else:
        return 
############################################################################################################################### 
# def magfetch(
#     # parameter = 'Bx',
#     start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
#     end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
#     magname = 'atu', 
#     is_verbose = False, 
#     tgopw = 'ResUseNoCom'
# ):
#     """
#     MAGFETCH 
#         Function to fetch data for a given magnetometer. Pulls from ai.cdas or DTU.

#         Arguments:
#             start, end   : datetimes of the start and end of sampled data range.
#             magname      : IAGA ID for magnetometer being sampled. e.g.: 'upn'
#             is_verbose   : Boolean for whether debugging text is printed.
#             tgopw        : Password for Tromsø Geophysical Observatory

#         Returns:
#             df           : pandas dataframe with columns ['UT', 'MAGNETIC_NORTH_-_H', 'MAGNETIC_EAST_-_E', 'VERTICAL_DOWN_-_Z']
#     """
#     # Pull password for TGO from local .txt file:
#     file = open("tgopw.txt", "r")
#     tgopw = file.read()
#     if(is_verbose): print('Found Tromsø Geophysical Observatory password: ' + tgopw)
#     if(magname in ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']): 
#         if(is_verbose): print('Collecting data for ' + magname + ' from TGO.')
#         # data = cdas.get_data(
#         #     'sp_phys',
#         #     'THG_L2_MAG_'+ magname.upper(),
#         #     start,
#         #     end,
#         #     ['thg_mag_'+ magname]
#         # )
#         foo = _download_remote_file(start, site=magname, verbose=is_verbose)
#         df = pd.read_csv('output/foo.txt', skiprows = 5, delim_whitespace=True)
#         # Convert the 'DD/MM/YYYY HH:MM:SS' column to datetime format
#         df['DD/MM/YYYY HH:MM:SS'] = df['DD/MM/YYYY'] + ' ' + df['HH:MM:SS']
#         df['UT'] = pd.to_datetime(df['DD/MM/YYYY HH:MM:SS'], format='%d/%m/%Y %H:%M:%S')

#         # Rename the columns
#         df.rename(columns={'X': 'MAGNETIC_NORTH_-_H', 'Y': 'MAGNETIC_EAST_-_E', 'Z': 'VERTICAL_DOWN_-_Z'}, inplace=True)

#         # Convert the dataframe to a dictionary matching ai.cdas output
#         data = {
#             'UT': df['UT'].to_numpy(),
#             'MAGNETIC_NORTH_-_H': df['MAGNETIC_NORTH_-_H'].to_numpy(),
#             'MAGNETIC_EAST_-_E': df['MAGNETIC_EAST_-_E'].to_numpy(),
#             'VERTICAL_DOWN_-_Z': df['VERTICAL_DOWN_-_Z'].to_numpy()
#         }
#     else:
#         data = cdas.get_data(
#             'sp_phys',
#             'THG_L2_MAG_'+ magname.upper(),
#             start,
#             end,
#             ['thg_mag_'+ magname]
#         )
#     if(is_verbose): print('Data for ' + magname + ' collected: ' + str(len(data['UT'])) + ' samples.')
#     return data

############################################################################################################################### 

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
    """
    MAGFIG 
        Function to create stacked plot of conjugate magnetometers.

        Arguments:
            parameter   : The parameter of interest - Bx, By, or Bz. North/South, East/West, and vertical, respectively.
            start, end   : datetimes of the start and end of plots
            maglist_a     : List of Arctic magnetometers. Default: ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']
            maglist_b     : Corresponding list of Antarctic magnetometers. Default: ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5']
            is_displayed   : Boolean for whether resulting figure is displayed inline. False by default.
            is_saved       : Boolean for whether resulting figure is saved to /output directory.
            events     : List of datetimes for events marked on figure. Empty by default.
            event_fontdict : Font dict for formatting of event labels. Default: {'size':20,'weight':'bold'}

        Returns:
            Figure of stacked plots for date in question, with events marked.
    """
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
            # data = cdas.get_data(
            #     'sp_phys',
            #     'THG_L2_MAG_'+ magname.upper(),
            #     start,
            #     end,
            #     ['thg_mag_'+ magname]
            # )
            # data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
            data = magfetch(start, end, magname)
                # Check the type of the data variable
            if isinstance(data, dict):
                # Extract the data for the specified parameter
                data = data[d[parameter]]

            x =data['UT']
            y =data[d[parameter]]
            y = reject_outliers(y) # Remove power cycling artifacts on, e.g., PG2.
            axs[idx].plot(x,y)#x, y)
            axs[idx].set(xlabel='Time', ylabel=magname.upper())

            if events is not None:
                # print('Plotting events...')
                trans      = mpl.transforms.blended_transform_factory(axs[idx].transData,axs[idx].transAxes)
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
                # data = cdas.get_data(
                #     'sp_phys',
                #     'THG_L2_MAG_'+ magname.upper(),
                #     start,
                #     end,
                #     ['thg_mag_'+ magname]
                # )
                # data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
                print('Pulling data.')
                data = magfetch(start, end, magname)
                if isinstance(data, dict):
                    # Extract the data for the specified parameter
                    data = data[d[parameter]]

                x =data['UT']
                y =data[d[parameter]]
                y = reject_outliers(y) # Remove power cycling artifacts on, e.g., PG2.
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
    fig.suptitle(str(start) + ' ' +  str(parameter), fontsize=30)   # Title the plot...
    if is_saved:
        print("Saving figure. " + fname)
        # fname = 'output/' +str(start) + '_' +  str(parameter) + '.png'
        fig.savefig(fname, dpi='figure', pad_inches=0.3)
    if is_displayed:
        return fig # TODO: Figure out how to suppress output here

###############################################################################################################################  
# Function to reject outliers. We'll need this to eliminate power cycling artifacts in the magnetometer plots.
def reject_outliers(y):   # y is the data in a 1D numpy array
    """
        Function to reject outliers from a 1D dataset.

        Arguments:
            y      : 1D numpy array

        Returns:
            array with outliers replaced with NaN
    """
    mean = np.mean(y)
    sd = np.std(y)
    final_list = np.copy(y)
    for n in range(len(y)):
        final_list[n] = y[n] if y[n] > mean - 3 * sd else np.nan
        final_list[n] = final_list[n] if final_list[n] < mean + 5 * sd else np.nan
    return final_list

############################################################################################################################### 
# # MAGSPECT Function to create power plots for conjugate magnetometers.

# def magspect(
#     parameter='Bx',
#     start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
#     end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
#     maglist_a=['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],
#     maglist_b=['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],
#     is_displayed=False,
#     is_saved=True,
#     is_verbose=False,
#     events=None,
#     event_fontdict={'size': 20, 'weight': 'bold'},
#     myFmt=mdates.DateFormatter('%H:%M')
# ):
#     """
#     Function to create power plots for conjugate magnetometers.

#     Arguments:
#         parameter: The parameter of interest - Bx, By, or Bz. North/South, East/West, and vertical, respectively.
#         start, end: datetimes of the start and end of plots
#         maglist_a: List of Arctic magnetometers. Default: ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']
#         maglist_b: Corresponding list of Antarctic magnetometers. Default: ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5']
#         is_displayed: Boolean for whether resulting figure is displayed inline. False by default.
#         is_saved: Boolean for whether resulting figure is saved to /output directory.
#         events: List of datetimes for events marked on figure. Empty by default.
#         event_fontdict: Font dict for formatting of event labels. Default: {'size': 20, 'weight': 'bold'}
#         myFmt: Date formatter. By default: mdates.DateFormatter('%H:%M')

#     Returns:
#         Figure of stacked plots for date in question, with events marked.
#     """

#     # Magnetometer parameter dict so that we don't have to type the full string:
#     d = {'Bx': 'MAGNETIC_NORTH_-_H', 'By': 'MAGNETIC_EAST_-_E', 'Bz': 'VERTICAL_DOWN_-_Z'}
#     if is_saved:
#         fname = 'output/' + str(start) + '_' + str(parameter) + '.png'
#         if os.path.exists(fname):
#             print('Looks like ' + fname + ' has already been generated.')
#             return

#         # raise Exception('This file has already been generated.')
#     fig, axs = plt.subplots(len(maglist_a), 2, figsize=(25, 25), constrained_layout=True)
#     print('Plotting data for ' + str(len(maglist_a)) + ' magnetometers: ' + str(start))

#     for maglist, side in zip([maglist_a, maglist_b], ['Arctic', 'Antarctic']):
#         for idx, magname in enumerate(maglist):
#             print('Plotting data for ' + side + ' magnetometer #' + str(idx+1) + ': ' + magname.upper())

#             try:
#                 data = magfetch(start, end, magname, is_verbose=is_verbose)
#                 # print(type(data['UT']))

#                 # data['UT'] = pd.to_datetime(data['UT'])  # Converting UT column to datetime object

#                 x = data['UT']
#                 y = data[d[parameter]]
#                 y = reject_outliers(y)
#                 # Interpolate NaNs away:
#                 df = pd.DataFrame(y, x)
#                 df = df.interpolate('linear')
#                 y = df[0].values

#                 xlim = [start, end]

#                 # axs[idx, side == 'Arctic'].xaxis.set_major_formatter(myFmt)
#                 # axs[idx, side == 'Arctic'].gca().xaxis.set_major_formatter(myFmt)
#                 # axs[:, side == 'Arctic'].gca().xaxis.set_major_formatter(myFmt)
#                 # axs[idx, side == 'Arctic'].set_major_formatter(myFmt)
#                 # axs[idx, side == 'Arctic'].gca().set_major_formatter(myFmt)

#                 f, t, Zxx = stft(y - np.mean(y), fs=1, nperseg=256, noverlap=128)
#                 dt_list = [start + dt.timedelta(seconds=ii) for ii in t]

#                 # cmap = axs[idx, side == 'Arctic'].pcolormesh(dt_list, f * 1000., np.abs(Zxx) * np)
#                 cmap=axs[idx, 0].pcolormesh(dt_list, f*1000., np.abs(Zxx)*np.abs(Zxx), vmin=0, vmax=0.5)
                

#                 axs[idx, 0].set(xlim=xlim)
#                 axs[idx, 0].set(ylim=[0,20])
#                 axs[idx, 0].xaxis.set_major_formatter(myFmt)
#                 axs[idx, 0].set_title('STFT Power Spectrum: ' + magname.upper())
#                 axs[idx, 0].set_ylabel('Frequency [mHz]')
#                 axs[idx, 0].set_xlabel('UT')
                
#                 axs[idx, 1].set(xlim=xlim)
#                 axs[idx, 1].set(ylim=[0,20])
#                 axs[idx, 1].xaxis.set_major_formatter(myFmt)
#                 axs[idx, 1].set_title('STFT Power Spectrum: ' + magname.upper())
#                 axs[idx, 1].set_ylabel('Frequency [mHz]')
#                 axs[idx, 1].set_xlabel('UT')
                
#                 if events is not None:
#                     # print('Plotting events...')
#                     trans      = mpl.transforms.blended_transform_factory(axs[idx,0].transData,axs[idx,0].transAxes)
#                     for event in events:
#                         evt_dtime   = event.get('datetime')
#                         evt_label   = event.get('label')
#                         evt_color   = event.get('color','0.4')

#                         axs[idx, 0].axvline(evt_dtime,lw=1,ls='--',color=evt_color)
#                         if evt_label is not None:
#                             axs[idx, 0].text(evt_dtime,0.01,evt_label,transform=trans,
#                                     rotation=90,fontdict=event_fontdict,color=evt_color,
#                                     va='bottom',ha='right')

#             except Exception as e:
#                 print(e)
#                 continue
#         fig.suptitle(str(start) + ' ' +  str(parameter), fontsize=30)   # Title the plot...
#         if is_saved:
#             fname = 'output/PowerSpectrum_' +str(start) + '_' +  str(parameter) + '.png'
#             print("Saving figure. " + fname)
#             fig.savefig(fname, dpi='figure', pad_inches=0.3)
#         if is_displayed:
#             return fig 


def magspect(
    parameter = 'Bx',
    start = datetime.datetime(2016,1,25,0,0), 
    end = datetime.datetime(2016,1,25,9,0), 
    maglist_a = ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],
    maglist_b = ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],
    is_displayed = False,
    is_saved = True, 
    is_verbose = False,
    events=None, event_fontdict = {'size':20,'weight':'bold'}, 
    myFmt = mdates.DateFormatter('%H:%M')
):
    """
       Function to create power plots for conjugate magnetometers.

        Arguments:
            parameter   : The parameter of interest - Bx, By, or Bz. North/South, East/West, and vertical, respectively.
            start, end  : datetimes of the start and end of plots
            maglist_a     : List of Arctic magnetometers. Default: ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']
            maglist_b     : Corresponding list of Antarctic magnetometers. Default: ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5']
            is_displayed   : Boolean for whether resulting figure is displayed inline. False by default.
            is_saved       : Boolean for whether resulting figure is saved to /output directory.
            events     : List of datetimes for events marked on figure. Empty by default.
            event_fontdict : Font dict for formatting of event labels. Default: {'size':20,'weight':'bold'}
            myFmt        : Date formatter. By default: mdates.DateFormatter('%H:%M')

        Returns:
            Figure of stacked plots for date in question, with events marked.
    """
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
            # data = cdas.get_data(
            #     'sp_phys',
            #     'THG_L2_MAG_'+ magname.upper(),
            #     start,
            #     end,
            #     ['thg_mag_'+ magname]
            # )
            # data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
            data = magfetch(start, end, magname, is_verbose = is_verbose)
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
                trans      = mpl.transforms.blended_transform_factory(axs[idx,0].transData,axs[idx,0].transAxes)
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
            # data = cdas.get_data(
            #     'sp_phys',
            #     'THG_L2_MAG_'+ magname.upper(),
            #     start,
            #     end,
            #     ['thg_mag_'+ magname]
            # )
            # data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
            print('Pulling data...')
            data = magfetch(start, end, magname, is_verbose = is_verbose)
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
                trans      = mpl.transforms.blended_transform_factory(axs[idx,1].transData,axs[idx,1].transAxes)
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
    fig.suptitle(str(start) + ' ' +  str(parameter), fontsize=30)   # Title the plot...
    if is_saved:
        fname = 'output/PowerSpectrum_' +str(start) + '_' +  str(parameter) + '.png'
        print("Saving figure. " + fname)
        fig.savefig(fname, dpi='figure', pad_inches=0.3)
    if is_displayed:
        return fig # TODO: Figure out how to suppress output here
        
        
############################################################################################################################### 
# # MAGDF Function to create multi-indexable dataframe of all mag parameters for a given period of time. 


def magdf(
    start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
    end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
    maglist_a = ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],  # Arctic magnetometers
    maglist_b = ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],  # Antarctic magnetometers
    is_saved = False, 
    ):
    """
       Function to create power plots for conjugate magnetometers.

        Arguments:
            start, end   : datetimes of the start and end of plots
            maglist_a     : List of Arctic magnetometers. Default: ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']
            maglist_b     : Corresponding list of Antarctic magnetometers. Default: ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5']
            is_saved       : Boolean for whether resulting dataframe is saved to /output directory.
            is_verbose    : Boolean for whether debugging text is printed. 

        Returns:
            Dataframe of Bx, By, Bz for each magnetometer in list.
    """
    # Magnetometer parameter dict so that we don't have to type the full string:
    d = {'Bx': 'MAGNETIC_NORTH_-_H', 'By': 'MAGNETIC_EAST_-_E', 'Bz': 'VERTICAL_DOWN_-_Z'}

    d_i = dict((v, k) for k, v in d.items()) # inverted mapping for col renaming later
    if is_saved:
        fname = 'output/' +str(start) + '_' + '.csv'
        if os.path.exists(fname):
            if(is_verbose): print('Looks like ' + fname + ' has already been generated.')
            return 
    UT = pd.date_range(start, end, freq ='S')   # preallocate time range
    full_df = pd.DataFrame(UT, columns=['UT'])   # preallocate dataframe
    full_df['UT'] = full_df['UT'].astype('datetime64[s]') # enforce 1s precision
    full_df['Magnetometer'] = ""
    for mags in [maglist_a, maglist_b]:
        for idx, magname in enumerate(mags):   # For each magnetometer, pull data and merge into full_df:
            if(is_verbose): print('Pulling data for magnetometer: ' + magname.upper())
            try:                
                # data = cdas.get_data(
                #     'sp_phys',
                #     'THG_L2_MAG_'+ magname.upper(),
                #     start,
                #     end,
                #     ['thg_mag_'+ magname]
                # )
                # data['UT'] = pd.to_datetime(data['UT'])# unit='s')
                # df = pd.DataFrame(data)
                df = magfetch(start, end, magname)
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
        if(is_verbose): print('Saving as a CSV.')
        full_df.to_csv(fname)
    # print(full_df)
    return full_df 
    
############################################################################################################################### 
# #  FILL_NAN: Function to eliminate NaN values from a 1D numpy array.

def fill_nan(y):
    """
        Fit a linear regression to the non-nan y values

        Arguments:
            y      : 1D numpy array with NaNs in it

        Returns:
            Same thing; no NaNs.
    """
    
    # Fit a linear regression to the non-nan y values

    # Create X matrix for linreg with an intercept and an index
    X = np.vstack((np.ones(len(y)), np.arange(len(y))))

    # Get the non-NaN values of X and y
    X_fit = X[:, ~np.isnan(y)]
    y_fit = y[~np.isnan(y)].reshape(-1, 1)

    # Estimate the coefficients of the linear regression
    beta = np.linalg.lstsq(X_fit.T, y_fit)[0]

    # Fill in all the nan values using the predicted coefficients
    y.flat[np.isnan(y)] = np.dot(X[:, np.isnan(y)].T, beta)
    return y

############################################################################################################################### 
# #  WAVEPWR: Function to determine Pc5 (by default) wave power for a given magnetometer, parameter and time frame.

def wavepwr(station_id, 
            parameter,         # Bx, By or Bz
            start, 
            end, 
            f_lower = 2.5,        # frequency threshold in mHz 
            f_upper = 3,     # frequency threshold in mHz
            is_verbose = False
           ):
    """
         Function to determine Pc5 (by default) wave power for a given magnetometer, parameter and time frame.

        Arguments: 
               station_id      : Station ID in lowercase, e.g., 'atu', 'pg4'
               parameter        : 'Bx', 'By' or 'Bz'
               start, end      : datetimes of interval
               f_lower, f_upper : Range of frequencies of interest in mHz.
               is_verbose      : Print details of calculation. False by default. 

        Returns:
               pwr        : Calculated wave power in range of interest. 
    """
    magname = station_id.lower()#'stf'
    d = {'Bx':'MAGNETIC_NORTH_-_H', 'By':'MAGNETIC_EAST_-_E','Bz':'VERTICAL_DOWN_-_Z'}
    # print(magname)
    try:
        if(is_verbose): print('Checking wave power for magnetometer ' + magname.upper() + ' between ' + str(start) + ' and ' + str(end) + '.')
        data = cdas.get_data(
            'sp_phys',
            'THG_L2_MAG_'+ magname.upper(),
            start,
            end,
            ['thg_mag_'+ magname]
        )
        # data = magfetch(start, end, magname)
        # data['UT'] = pd.to_datetime(data['UT'])#, unit='s')
        x =data['UT']
        y =data[d[parameter]]


        y = reject_outliers(y) # Remove power cycling artifacts on, e.g., PG2.
        y = fill_nan(y)
        y = y - np.nanmean(y)  # Detrend
        # print(y)

        dt = (x[1] - x[0]).seconds
        fs = 1 / dt

        datos = y

        # nblock = 1024
        # overlap = 128
        nblock = 60
        overlap = 30
        win = hann(nblock, True)

        # f, Pxxf = welch(datos, fs, window=win, noverlap=overlap, nfft=nblock, return_onesided=True, detrend=False)
        f, Pxxf = welch(datos, fs, window=win, return_onesided=True, detrend=False)

#       plt.semilogy(f, Pxxf, '-o')
#       plt.xlabel('Frequency [Hz]')
#       plt.ylabel('PSD [V**2/Hz]')


#       plt.grid()
        # plt.show()
        # pwr = np.trapz(Pxxf[((f>=f_lower/1000) & (f_upper<=3/1000))])
        # print(f)
        # print(Pxxf)
        # pwr = sum(Pxxf[((f>=f_lower/1000) & (f_upper<=3/1000))])
        pwr = Pxxf[3]
        if(is_verbose): print(Pxxf[((f>=f_lower/1000) & (f_upper<=3/1000))])
        if(is_verbose): print(magname.upper() + ': The estimated power from ' + str(f_lower) + ' mHz to '+ str(f_upper) + ' mHz is ' + str(pwr) + ' nT/Hz^(1/2)')
        # print(pwr)
        return pwr
    except Exception as e:
        print(e)
        # print('Window length: ' + str(len(win)) +'; '+'Signal length: ' + str(len(y))) # usually this is the issue.
        return 'Error'
    

############################################################################################################################### 

def wavefig(
    stations, # dataframe 
    parameter = 'Bx',
    start = datetime.datetime(2016, 1, 24, 0, 0, 0), 
    end = datetime.datetime(2016, 1, 25, 0, 0, 0), 
    maglist_a = ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb'],
    maglist_b = ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5'],
    is_displayed = True,
    is_saved = False, 
    is_verbose = False
):
    """
    WAVEFIG 
        Function to create wave power plot for a given set of magnetometers. 

        Arguments:
            stations    : Dataframe of stations with columns IAGA, AACGMLAT, AACGMLON
            parameter   : The parameter of interest - Bx, By, or Bz. North/South, East/West, and vertical, respectively.
            start, end   : datetimes of the start and end of plots
            maglist_a     : List of Arctic magnetometers. Default: ['upn', 'umq', 'gdh', 'atu', 'skt', 'ghb']
            maglist_b     : Corresponding list of Antarctic magnetometers. Default: ['pg0', 'pg1', 'pg2', 'pg3', 'pg4', 'pg5']
            is_displayed   : Boolean for whether resulting figure is displayed inline. False by default.
            is_saved       : Boolean for whether resulting figure is saved to /output directory.
            is_verbose     : Boolean for whether debugging text is printed.

        Returns:
            Figure of stacked plots for date in question, with events marked.
    """
    
    foo = stations.copy()
    foo['WAVEPWR'] = foo.apply(lambda row : wavepwr(row['IAGA'], parameter = parameter, start = start, end = end, is_verbose = is_verbose), axis = 1)#wavepwr(foo['IAGA'])
    foo['HEMISPHERE'] = np.sign(foo.AACGMLAT)
    foo.HEMISPHERE = foo['HEMISPHERE'].map({1:'Arctic', -1:'Antarctic', 0:'Error'})
    foo['ABSLAT'] = abs(foo.AACGMLAT)
    px.scatter(foo.WAVEPWR, abs(foo.AACGMLAT))
    foo = foo.sort_values('ABSLAT')
    
    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    x = foo.query("HEMISPHERE == 'Antarctic'")['ABSLAT'].to_list()
    y = foo.query("HEMISPHERE == 'Antarctic'")['WAVEPWR'].to_list()
    labels = foo.query("HEMISPHERE == 'Antarctic'")['IAGA'].to_list()

    # Add traces
    fig.add_trace(
        go.Scatter(x= x,y = y, text = labels, name="Antarctic"),
        secondary_y=False,
    )

    x = foo.query("HEMISPHERE == 'Arctic'")['ABSLAT'].to_list()
    y = foo.query("HEMISPHERE == 'Arctic'")['WAVEPWR'].to_list()
    labels = foo.query("HEMISPHERE == 'Arctic'")['IAGA'].to_list()

    fig.add_trace(
        go.Scatter(x=x, y=y, text = labels, name="Arctic", 
                   # line_shape = 'spline'
                  ),
        secondary_y=True,
    )

    # Add figure title
    fig.update_layout(
        title_text=str(parameter) + " Wave Power: " + str(start) + " to " + str(end)
    )

    # Set x-axis title
    fig.update_xaxes(title_text="Latitude")

    # Set y-axes titles
    fig.update_yaxes(title_text="<b>antarctic</b> wave power", secondary_y=False)
    fig.update_yaxes(title_text="<b>arctic</b> wave power", secondary_y=True)
    
    if(is_saved): 
        if(is_verbose): print("Saving figure.")
        fname = 'output/WavePower_' +str(start) + '_to_' + str(end) + '_' +  str(parameter) + '.png'
        print("Saving figure. " + fname)
        # fig.savefig(fname, dpi='figure', pad_inches=0.3)
        fig.write_image(fname)
    if is_displayed:
        # fig.show()
        return fig # TODO: Figure out how to suppress output here
############################################################################################################################### 

def _download_remote_file(datetime: dt.datetime, site='tro', verbose=False):
    """Downloads file for a given date/site/orientation to the configured storage path

    Args:
        datetime (dt.datetime): date/time to download
        site (str, optional): Which site code
        orientation (str, optional): sensor orientation ('XYZ', 'HEZ')
    """
    def _build_request(datetime, site='tro'):
        request_dict = {#'site':config.site_codes[site],
                        'year':datetime.strftime('%Y'),
                        'month':'{}'.format(datetime.month),
                        'day':'{}'.format(datetime.day),
                        'res':'10sec',
                        'pwd': 'ResUseNoCom',
                        'format':'XYZhtml',
                        'comps':'SXYZ',
                        'getdata':' Get Data '}
        return request_dict
    
    # initiate request on remote server
    input_url = 'http://flux.phys.uit.no/cgi-bin/mkascii.cgi?'

    request_dict = _build_request(datetime, site=site)
    init_resp = requests.get(input_url, params=request_dict, timeout=(30,300))
    
    # check if header and footer are as expected
    if init_resp.text.startswith('<pre>\n') and init_resp.text.endswith('</pre>\n'):
        savedata = init_resp.text.lstrip('<pre>\n').rstrip('</pre>\n') + '\n'
        # find a place to put the data
        savefilepath = 'output/' + datetime.strftime(datetime) + '_' + site + '_' + sample_rate +'.txt'
        if(verbose): print(savefilepath)
        try:
            with open(savefilepath, 'w') as savefile:
                savefile.write(savedata)
        except FileNotFoundError:
            # os.makedirs(datetime.strftime(config.fpath_format.format(base=config.datapath_local)))
            with open(savefilepath, 'w') as savefile:
                savefile.write(savedata)
        return savedata
    else:
        return 