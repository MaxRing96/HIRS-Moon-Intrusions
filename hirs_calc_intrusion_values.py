import warnings
# ignore some unimportant warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', 
message='The unit of the quantity is stripped when downcasting to ndarray.')

import numpy as np
import pandas as pd
import hirs_src as src
import sys

print('\n**********')

# read in the filepath (input-filename) of the file with the 
# suspected intrusion given by the argument in the command line
intrusion_file = sys.argv[1]

# retrieve date of the intrusion out of the input-filename 
# string to use it in the output-filename
date_of_intrusion = (
    f'{intrusion_file[40:44]}-'
    f'{intrusion_file[45:47]}-'
    f'{intrusion_file[48:50]}'
)
print('Date of intrusion: ',date_of_intrusion)

# retrieve satellite name from input-filename 
SATELLITE = str(intrusion_file[28:34])
print(f'\nSatellite: {SATELLITE}')

# set up hirs channels 1-19
CHANNELS = np.arange(1,20)

# set up hirs data handler
read_hirs = src.get_hirs_reader(SATELLITE)

# read in hirs-file of moon intrusion as xarray dataset
(lines, extra) = read_hirs.read(intrusion_file)
ds = read_hirs.as_xarray_dataset(lines)

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
    
# set up empty lists for data storage
ch_central_waves = []
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

    # find index and timestamp of scanline with minimum mean dsv counts
    # (= scnanline of moon intrusion)
    intrusion_timeindex = counts_scanlines.mean(dim="scanpos").argmin()
    intrusion_timestamp = str(counts_scanlines.isel(
        time=intrusion_timeindex).time.values)[11:19]

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
        # If no meta file is found, raise error message but continue with loop.
    try:
        sat_meta_data = pd.read_csv(f'meta_files/meta_{SATELLITE}.txt', comment='#').set_index('channel')

        ch_wavenumber = sat_meta_data.loc[CHANNELS[j]]['wavenumber']
        ch_central_waves.append(ch_wavenumber)
        corr_fac1 = sat_meta_data.loc[CHANNELS[j]]['correction_factor1']
        corr_fac2 = sat_meta_data.loc[CHANNELS[j]]['correction_factor2']

        bb_radiance = src.calc_bb_flux(ch_wavenumber*10**2,\
        [corr_fac1,corr_fac2],bb_temp_mean)

        bb_rad_all.append(round(bb_radiance,2))

    except FileNotFoundError:
        print(f'\n*** No meta data file for {SATELLITE} found ***\n'
        '*** --> Blackbody radiance could not be calculated ***')
        continue

# get 3D position of the satellite
longitude, latitude, altitude = src.get_position(ds,intrusion_timeindex)

### Save data to csv file ###

# build dictionary with calculated data
data = {
    'channel': CHANNELS,
    'central_wavenumber' : ch_central_waves,
    'intrusion_scan_range': intrusion_range_all,
    'bb_temp_mean': round(bb_temp_mean,4),
    'bb_temp_std': round(bb_temp_std,4),
    'moon_counts_mean': moon_counts_mean_all,
    'moon_counts_std': moon_counts_std_all,
    'dsv_counts_mean': dsv_counts_mean_all,
    'dsv_counts_std': dsv_counts_std_all,
    'bb_counts_mean': bb_counts_mean_all,
    'bb_counts_std': bb_counts_std_all,
    'longitude': longitude.tolist(),
    'latitude': latitude.tolist(),
    'altitude': altitude.tolist(),
}

# convert dictionary to panda dataframe
df = pd.DataFrame(data)

# if blackbody radiance could be calculated  
# insert the data into the dataframe at column 10 
if len(bb_rad_all) != 0:
    df.insert(10, "bb_rad", bb_rad_all, True)

# set path to which the output csv file should be saved to
# default path: local dircetory
PATH = ''
FILENAME = f'hirs_moon_intrusion_{SATELLITE}_{date_of_intrusion}_{intrusion_timestamp}_calculations.csv'

# convert dataframe to csv and save it 
df = df.set_index('channel')
df.to_csv(PATH+FILENAME)

print('')
print('Calculation results:')
print('')
print(df)
print('')
print(f'Output saved to:    {PATH+FILENAME}')
print('****************')
print('')
print("Now go to Horizon webpage https://ssd.jpl.nasa.gov/horizons.cgi and get the angular diameter and the phase angle of the moon.")
print('')
print("You will need the following input variables:")
print("Ephemeris Type: OBSERVER")
print("Target Body: Moon[Luna][301]")
print("Observer Location: Topocentric (lon, lat, alt): ", (df.longitude[1], df.latitude[1], df.altitude[1]))
print('Time span: Date of intrusion "',date_of_intrusion, intrusion_timestamp, '", "Date_of_intrusion+1min", Step=1m')
print("Table Settings: QUANTITIES= 13 (Target angular diamenter), 24 (Sun-Target-Observer ~PHASE angle)")
print("Display/Output: default")
print('****************')
print("")
#EOF