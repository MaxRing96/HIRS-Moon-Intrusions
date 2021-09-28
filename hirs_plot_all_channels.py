import warnings
# ignore some unimportant warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', 
message='The unit of the quantity is stripped when downcasting to ndarray.')

import numpy as np
import pandas as pd
from typhon.plots import styles
import matplotlib.pyplot as plt
import sys
import hirs_src as src

styles.use('typhon')

print('\n**********')

# read in the filepath (input-filename) of the file with the 
# suspected intrusion given by the argument in the command line
intrusion_file = sys.argv[1]

# read in the file with the channel-specified moon intrusion scan ranges.
# Raise error if file not found and exit.
try:
    intr_scan_ranges = pd.read_csv('config_intr_scan_ranges.txt').set_index('channel')
except FileNotFoundError:
    print('')
    print('*** Textfile config_intr_scan_ranges.txt with the specified scan ranges of '
    'the moon intrusion for each channel could not be found! ***')
    print('')
    raise

# set up data handler for hirs 
SATELLITE = str(intrusion_file[28:34])
print('\nSatellite: ',SATELLITE)
read_hirs = src.get_hirs_reader(SATELLITE)


# read in hirs file of moon intrusion as xarray dataset
(lines, extra) = read_hirs.read(intrusion_file)
ds = read_hirs.as_xarray_dataset(lines)

# set all the text in plot to size 10 
plt.rcParams['font.size'] = '10'

# plotting loop:
# plots a 5x4 matrix of plots each containing the dsv counts
# for the individual channel during the intrusion.
(fig, axs) = plt.subplots(nrows=5,ncols=4,sharex=True,figsize=(14,12))

# this loop plots channels 1-5
for i in range(0,5):
    
    channels = np.arange(1,6)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    # find scanline of intrusion by searching for the 
    # scanline with the minimum average of counts:
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)

    
    # retrieve the scan time of the intrusion to use it in the 
    # plot title and filename
    scan_time = str(scanline_of_intrusion_counts.time.values)[11:19]

    # plot dsv counts measured during intrusion by channel i:
    axs[i,0].grid(axis='x',color='grey',linestyle='--',linewidth=0.4)
    axs[i,0].plot(scanline_of_intrusion_counts.scanpos,
                scanline_of_intrusion_counts)
    axs[i,0].set_title(f'channel {channels[i]}\n'
                       f'scan time: {scan_time}',
                       fontsize=8)
    axs[i,0].set_xlim(10,56)

    # plot the lower and upper bound of the specified
    # scan-range (specified in the 'config_intr_scan_ranges.txt' file) 
    # of the intrusion as vertical lines
    axs[i,0].vlines(
        x=intr_scan_ranges.loc[channels[i]]['min_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        )

    axs[i,0].vlines(
        x=intr_scan_ranges.loc[channels[i]]['max_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        label='intrusion range',
        )

# this loop plots channels 6-10
for i in range(0,5):
    
    channels = np.arange(6,11)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    axs[i,1].grid(axis='x',color='grey',linestyle='--',linewidth=0.4)
    axs[i,1].plot(scanline_of_intrusion_counts.scanpos,
                scanline_of_intrusion_counts)
    axs[i,1].set_title(f'channel {channels[i]}\n'
                       f'scan time: {str(scanline_of_intrusion_counts.time.values)[11:19]}',
                       fontsize=8)
    axs[i,1].set_xlim(10,56)

    axs[i,1].vlines(
        x=intr_scan_ranges.loc[channels[i]]['min_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        )

    axs[i,1].vlines(
        x=intr_scan_ranges.loc[channels[i]]['max_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        label='intrusion range',
        )
    
# this loop plots channels 11-15
for i in range(0,5):
    
    channels = np.arange(11,16)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    axs[i,2].grid(axis='x',color='grey',linestyle='--',linewidth=0.4)
    axs[i,2].plot(scanline_of_intrusion_counts.scanpos,
                scanline_of_intrusion_counts)
    axs[i,2].set_title(f'channel {channels[i]}\n'
                       f'scan time: {str(scanline_of_intrusion_counts.time.values)[11:19]}',
                       fontsize=8)
    axs[i,2].set_xlim(10,56)

    axs[i,2].vlines(
        x=intr_scan_ranges.loc[channels[i]]['min_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        )

    axs[i,2].vlines(
        x=intr_scan_ranges.loc[channels[i]]['max_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        label='intrusion range',
        )

# this loop plots channels 16-19
# (last plot is empty) 
for i in range(0,4):

    channels = np.arange(16,20)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    axs[i,3].grid(axis='x',color='grey',linestyle='--',linewidth=0.4)
    axs[i,3].plot(scanline_of_intrusion_counts.scanpos,
                scanline_of_intrusion_counts)
    axs[i,3].set_title(f'channel {channels[i]}\n'
                       f'scan time: {str(scanline_of_intrusion_counts.time.values)[11:19]}',
                       fontsize=8)
    axs[i,3].set_xlim(10,56)

    axs[i,3].vlines(
        x=intr_scan_ranges.loc[channels[i]]['min_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        )

    axs[i,3].vlines(
        x=intr_scan_ranges.loc[channels[i]]['max_scanpos'],
        ymin=scanline_of_intrusion_counts.min(),
        ymax=scanline_of_intrusion_counts.max(),
        color='red',
        linewidth=1,
        linestyle='--',
        label='intrusion range',
        )
    
    handles, labels = axs[i,3].get_legend_handles_labels()

# set general x- and y-labels for matrix plot as well as the label
fig.text(0.5, -0.02, 'scanpos', ha='center',fontsize=13,fontweight="bold")
fig.text(-0.02, 0.5, 'counts', va='center', rotation='vertical',fontsize=13,fontweight="bold")
fig.legend(handles, labels, loc='lower right', prop={'size': 12})

# retrieve the date of the intrusion out of the timestamp 
# in the data to use it in the plot title and output-filename
date_of_intrusion = str(scanline_of_intrusion_counts.time.values)[:10]
print('Date of intrusion: ',date_of_intrusion)

# set plot title
fig.suptitle(f"POSSIBLE MOON INTRUSION for {SATELLITE}, DATE: {date_of_intrusion}\n"
             f"granule file: {intrusion_file}",fontsize=12)

# adjust spacing in between the plots   
fig.tight_layout(pad=1.2,w_pad=0)

# set path to which the plot should be saved to
# default path: local dircetory
PATH = 'hirs_moon_intrusion_suspects_plots/'
FILE = f'hirs_moon_intrusion_{SATELLITE}_{date_of_intrusion}_{scan_time}_allchannels_dsv.pdf'
plt.savefig(PATH+FILE,bbox_inches='tight')


# set up hirs channels 1-19
CHANNELS = np.arange(1,20)
# loop through all channels
for j in range(0,len(CHANNELS)):

    # read in counts of all deep-space-view (dsv) & 
    # blackbody-target (bb) scanlines for channel j
    counts_scanlines = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),
        channel=CHANNELS[j])
    counts_scanlines_bb = ds[{"time":ds["scantype"]==3}]["counts"].sel(
        scanpos=range(10,56),
        channel=CHANNELS[j])

    # find index and timestamp of scanline with minimum mean dsv counts
    # (= scnanline of moon intrusion)
    intrusion_timeindex = counts_scanlines.mean(dim="scanpos").argmin()
    intrusion_timestamp = str(counts_scanlines.isel(
        time=intrusion_timeindex).time.values)[11:19]

# get 3D position of the satellite
longitude, latitude, altitude = src.get_position(ds,intrusion_timeindex)

print('\n**************')
print('Plot saved to: ',PATH+FILE)
print('')
print('You can adjust the ranges of the suspected moon intrusion in the file "config_intr_scan_ranges.txt".')
print('')
print('****************')
print('')
print("Next step is to go to Horizon webpage https://ssd.jpl.nasa.gov/horizons.cgi and get the ANGULAR DIAMETER and the PHASE ANGLE of the moon.")
print('')
print("You will need the following input variables:")
print("Ephemeris Type: OBSERVER")
print("Target Body: Moon[Luna][301]")
print("Observer Location: Topocentric (lon, lat, alt): ", longitude, latitude, altitude)
print('Time span: Date of intrusion "',date_of_intrusion, intrusion_timestamp, '", "Date_of_intrusion+1min", Step=1m')
print("Table Settings: QUANTITIES= 13 (Target angular diamenter), 24 (Sun-Target-Observer ~PHASE angle)")
print("Display/Output: default")
print('')
print('****************')
print('Then proceed with the calculations in hirs_calc_intrusion_values.py')
print("")
#EOF