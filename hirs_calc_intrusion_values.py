import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import hirs_src as src
import sys

print('**********')

# read in the filepath of the file with the suspected intrusion
# given by the argument in the command line
intrusion_file = sys.argv[1]

# read satellite name from file path
SATELLITE = str(intrusion_file[28:34])
print(f'\nSatellite: {SATELLITE}')

# set up hirs channels 1-19
CHANNELS = np.arange(1,20)

# set up hirs data handler
read_hirs = src.get_hirs_reader(SATELLITE)

# read in hirs file of moon intrusion
(lines, extra) = read_hirs.read(intrusion_file)
ds = read_hirs.as_xarray_dataset(lines)

# read in the file with the channel-specified moon intrusion scan ranges
try:
    intr_scan_ranges = pd.read_csv('config_intr_scan_ranges.txt').set_index('channel')
except FileNotFoundError:
    print('')
    print('*** Textfile config_intr_scan_ranges.txt with the specified scan ranges of '
    'the moon intrusion for each channel could not be found! ***')
    print('')
    raise
    
# set up empty lists for data storage
intrusion_range_all = []
dsv_counts_mean_all = []
dsv_counts_std_all = []
moon_counts_mean_all = []
moon_counts_std_all = []
bb_counts_mean_all = []
bb_counts_std_all = []
bb_rad_all = []

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

    # find index of scanline with minimum mean counts
    # (= scnanline of moon intrusion)
    intrusion_timeindex = counts_scanlines.mean(dim="scanpos").argmin()

    # read in the prescriped range of the moon intrusion
    intrusion_range = [
        intr_scan_ranges.loc[CHANNELS[j]]['min_scanpos'],
        intr_scan_ranges.loc[CHANNELS[j]]['max_scanpos']
        ]

    ### calculation of counts ###
    dsv_counts_mean, dsv_counts_std = src.calc_mean_counts_dsv(
        counts_scanlines,intrusion_timeindex)
    moon_counts_mean, moon_counts_std = src.calc_mean_counts(
        counts_scanlines,intrusion_timeindex,intrusion_range)
    bb_counts_mean, bb_counts_std = src.calc_mean_counts(
        counts_scanlines_bb,intrusion_timeindex)

    #### calculation of blackbody temperature ###
    # number of blackbody temperature sensors
    nr_bb_sensors = len(ds.prt_number_iwt.values)

    bb_temp_mean, bb_temp_std = src.calc_bb_temp(ds,intrusion_timeindex,\
        nr_bb_sensors)   

    # append calculated values to lists
    intrusion_range_all.append(intrusion_range)
    dsv_counts_mean_all.append(round(dsv_counts_mean,2))
    dsv_counts_std_all.append(round(dsv_counts_std,2))
    moon_counts_mean_all.append(round(moon_counts_mean,2))
    moon_counts_std_all.append(round(moon_counts_std,2))
    bb_counts_mean_all.append(round(float(bb_counts_mean),2))
    bb_counts_std_all.append(round(bb_counts_std,2))

    ### calculation of blackbody radiance ###
        # try to read in txt file with meta information for 
        # the satellite and channel, to calculate the blackbody radiance.
        # If no meta file is found, raise error message. 
    try:
        sat_meta_data = pd.read_csv(f'meta_{SATELLITE}.txt').set_index('channel')

        ch_wavenumber = sat_meta_data.loc[CHANNELS[j]]['wavenumber']
        corr_fac1 = sat_meta_data.loc[CHANNELS[j]]['correction_factor1']
        corr_fac2 = sat_meta_data.loc[CHANNELS[j]]['correction_factor2']

        bb_radiance = src.calc_bb_flux(ch_wavenumber*10**2,\
        [corr_fac1,corr_fac2],bb_temp_mean)

        bb_rad_all.append(bb_radiance)

    except FileNotFoundError:
        print(f'\n*** No meta data file for {SATELLITE} found ***\n'
        '*** --> Blackbody radiance could not be calculated ***')
        continue

# get platform information about longitude and altitude
longitude, latitude, altitude = src.get_position(ds,intrusion_timeindex)

### save output to csv file ###
date_of_intrusion = (
    f'{intrusion_file[40:44]}-'
    f'{intrusion_file[45:47]}-'
    f'{intrusion_file[48:50]}'
)

PATH = 'moon_intrusion_calculations/'
FILENAME = f'hirs_moon_intrusion_{SATELLITE}_{date_of_intrusion}_calculations.csv'

data = {
    'channel': CHANNELS,
    'intrusion scan-range': intrusion_range_all,
    'bb temp mean': round(bb_temp_mean,4),
    'bb temp std': round(bb_temp_std,4),
    'moon counts mean': moon_counts_mean_all,
    'moon counts std': moon_counts_std_all,
    'dsv counts mean': dsv_counts_mean_all,
    'dsv counts std': dsv_counts_std_all,
    'bb counts mean': bb_counts_mean_all,
    'bb counts std': bb_counts_std_all,
    'bb rad': bb_rad_all,
    'longitude': longitude.tolist(),
    'latitude': latitude.tolist(),
    'altitude': altitude.tolist(),
}

df = pd.DataFrame(data).set_index('channel')
df.to_csv(PATH+FILENAME)

print('')
print('Calculation results:')
print('')
print(df)
print('')
print('****************')
print('All calculations finished!')
print(f'Output saved to:    {PATH+FILENAME}')

#EOF