import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import typhon
import typhon.datasets.tovs
from typhon.plots import styles
import matplotlib.pyplot as plt
import sys
import hirs_src as src

styles.use('typhon')

intrusion_file = sys.argv[1]

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
read_hirs = src.get_hirs_reader(SATELLITE)

# read in hirs file of moon intrusion
(lines, extra) = read_hirs.read(intrusion_file)
ds = read_hirs.as_xarray_dataset(lines)

plt.rcParams['font.size'] = '10'

# plotting
(fig, axs) = plt.subplots(nrows=5,ncols=4,sharex=True,figsize=(14,12))

for i in range(0,5):
    
    channels = np.arange(1,6)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    # find scanline of intrusion by searching for the 
    # scanline with the minimum average of counts:
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    # plot counts measured during intrusion by channel i:
    axs[i,0].grid(axis='x',color='grey',linestyle='--',linewidth=0.4)
    axs[i,0].plot(scanline_of_intrusion_counts.scanpos,
                scanline_of_intrusion_counts)
    axs[i,0].set_title(f'channel {channels[i]}\n'
                       f'scan time: {str(scanline_of_intrusion_counts.time.values)[11:19]}',
                       fontsize=8)
    axs[i,0].set_xlim(10,56)

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
    
for i in range(0,5):
    
    channels = np.arange(6,11)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    # find scanline of intrusion by searching for the 
    # scanline with the minimum average of counts:
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    # plot counts measured during intrusion by channel i:
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
    
for i in range(0,5):
    
    channels = np.arange(11,16)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    # find scanline of intrusion by searching for the 
    # scanline with the minimum average of counts:
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    # plot counts measured during intrusion by channel i:
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
    
for i in range(0,4):
    
    channels = np.arange(16,20)
    
    counts = ds[{"time":ds["scantype"]==1}]["counts"].sel(
        scanpos=range(10,56),channel=channels[i])
    
    # find scanline of intrusion by searching for the 
    # scanline with the minimum average of counts:
    timeindex_of_intrusion = counts.mean(dim="scanpos").argmin()
    scanline_of_intrusion_counts = counts.isel(time=timeindex_of_intrusion)
    
    # plot counts measured during intrusion by channel i:
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
    
fig.text(0.5, -0.02, 'scanpos', ha='center',fontsize=13,fontweight="bold")
fig.text(-0.02, 0.5, 'counts', va='center', rotation='vertical',fontsize=13,fontweight="bold")

fig.legend(handles, labels, loc='lower right', prop={'size': 12})

#axs[0,0].vlines(x=20,ymin=scanline_of_intrusion_counts.min(),ymax=scanline_of_intrusion_counts.max())

DATE_OF_INTRUSION = str(scanline_of_intrusion_counts.time.values)[:10]

fig.suptitle(f"POSSIBLE MOON INTRUSION for {SATELLITE}, DATE: {DATE_OF_INTRUSION}\n"
             f"granule file: {intrusion_file}",fontsize=12)

    
fig.tight_layout(pad=1.2,w_pad=0)

PATH = f'moon_intrusion_plots/{SATELLITE}/possible_intrusions/'
FILE = f'hirs_moon_intrusion_{SATELLITE}_{DATE_OF_INTRUSION}_allchannels_dsv.pdf'
plt.savefig(PATH+FILE,bbox_inches='tight')

print('\n**************')
print('Plot saved to: ',PATH+FILE)
