import warnings
# ignore some unimportant warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', 
message='The unit of the quantity is stripped when downcasting to ndarray.')

import numpy as np
import pandas as pd
import hirs_src as src
import sys
from typhon.physics import planck_wavenumber, radiance2planckTb

print('\n**********')

# read in the filepath (input-filename) of the file with the 
# suspected intrusion given by the argument in the command line
intrusion_file = sys.argv[1]

# read in angular diam[arcsec] given by the second argument in the command line
ang_diam = sys.argv[2]
#conversion from arcsec to deg
ang_diam = float(ang_diam)/3600

# read in phase angle[deg] given by the third argument in the command line
phase_angle = sys.argv[3]

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
read_hirs, instrument = src.get_hirs_reader(SATELLITE)

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
ch_wavelength = []
ch_fov = []
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
        lamb = float(ch_wavenumber)*10**(-4)
        ch_wavelength.append(np.divide(1,lamb))
        corr_fac1 = sat_meta_data.loc[CHANNELS[j]]['correction_factor1']
        corr_fac2 = sat_meta_data.loc[CHANNELS[j]]['correction_factor2']
        fov = sat_meta_data.loc[CHANNELS[j]]['fov']
        ch_fov.append(fov)

        bb_radiance = src.calc_bb_flux(ch_wavenumber*10**2,\
        [corr_fac1,corr_fac2],bb_temp_mean)

        bb_rad_all.append(round(bb_radiance,2))

    except FileNotFoundError:
        print(f'\n*** No meta data file for {SATELLITE} found ***\n'
        '*** --> Blackbody radiance could not be calculated ***')
        continue

# get 3D position of the satellite
longitude, latitude, altitude = src.get_position(ds,intrusion_timeindex)

# build dictionary with calculated data
data = {
    'timestamp' : str(date_of_intrusion)+ " "+str(intrusion_timestamp),
    'satellite' : SATELLITE,
    'instrument': instrument,
    'channel': CHANNELS,
    'central_wavenumber[1/cm]' : ch_central_waves,
    'wavelength[micron]' : ch_wavelength ,
    'fov[deg]': ch_fov,
    'longitude[deg]': longitude.tolist(),
    'latitude[deg]': latitude.tolist(),
    'altitude[km]': altitude.tolist(),
    'intrusion_scan_range': intrusion_range_all,
    'bb_temp_mean[K]': round(bb_temp_mean,4),
    'bb_temp_std[K]': round(bb_temp_std,4),
    'moon_counts_mean': moon_counts_mean_all,
    'moon_counts_std': moon_counts_std_all,
    'dsv_counts_mean': dsv_counts_mean_all,
    'dsv_counts_std': dsv_counts_std_all,
    'bb_counts_mean': bb_counts_mean_all,
    'bb_counts_std': bb_counts_std_all,

}

# convert dictionary to panda dataframe
df = pd.DataFrame(data)

# if blackbody radiance could be calculated  
# insert the data into the dataframe at column 10 
if len(bb_rad_all) != 0:
    df.insert(10, "bb_radiance[MJy/sr]", bb_rad_all, True)

#add ang-diam[deg] and phase-angle[deg] to df
df["ang_diam[deg]"]=ang_diam
df["phase_angle[deg]"]= phase_angle

#calc slope
df["slope"] = np.divide(df["bb_radiance[MJy/sr]"],(df["bb_counts_mean"]-df["dsv_counts_mean"]))
df["slope_sigma"] = np.divide(df["slope"]*np.sqrt(df["bb_counts_std"]**2+df["dsv_counts_std"]**2),(df["bb_counts_mean"]-df["dsv_counts_mean"]))

#calc radiance
df["radiance[MJy/sr]"]=np.divide((df["bb_radiance[MJy/sr]"]+df["slope"]*(df["moon_counts_mean"]-df["bb_counts_mean"])),(df["ang_diam[deg]"]**2))*(df["fov[deg]"]**2)/0.97
df["radiance_sigma[%]"]=abs(100*np.sqrt((np.divide(df["slope_sigma"],df["slope"]))**2+np.divide(((df["bb_counts_std"])**2+(df["moon_counts_std"])**2),((df["bb_counts_mean"]-df["moon_counts_mean"])**2)))*np.divide((df["moon_counts_mean"]-df["bb_counts_mean"]),(df["dsv_counts_mean"]-df["bb_counts_mean"])))

#speed of light [m/s]
co = 2.9979245e+08

#calc TB
df["brightness_temp[K]"] = radiance2planckTb((co*df["central_wavenumber[1/cm]"]*1E02),(df["radiance[MJy/sr]"]*1E-20))

# set path to which the output csv file should be saved to
# default path: local dircetory
PATH = 'moon_intrusion_calculations_output/'
FILENAME = f'hirs_moon_intrusion_{SATELLITE}_{date_of_intrusion}_{intrusion_timestamp}_calculations.csv'

# convert dataframe to csv and save it 
df.set_index("channel", inplace=True)
df.to_csv(PATH+FILENAME)

Tb_mean_lw = df.loc[0:6,"brightness_temp[K]"].mean()
Tb_std_lw = df.loc[0:6,"brightness_temp[K]"].std()
Tb_mean_sw = df.loc[12:15,"brightness_temp[K]"].mean()
Tb_std_sw = df.loc[12:15,"brightness_temp[K]"].std()
print('')
print('Calculation results:')
print('')
print(np.round(df["brightness_temp[K]"],2))
print('')
print("Mean brightness Temperature LW Channel 2-7: ", np.round(Tb_mean_lw,2),"+/-", np.round(Tb_std_lw,2))
print("Mean brightness Temperature SW Channel 13-16: ", np.round(Tb_mean_sw,2),"+/-", np.round(Tb_std_sw,2))
print("Phase angle: ", df.loc[1,"phase_angle[deg]"])
print(f'Output saved to:    {PATH+FILENAME}')
print('****************')
print("")
#EOF