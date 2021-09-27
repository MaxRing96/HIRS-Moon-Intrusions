import numpy as np
import typhon
import typhon.datasets.tovs
from typhon.physics import planck_wavenumber, radiance2planckTb

# set up the fitting hirs data handler (for HIRS2, HIRS3 or HIRS4)
# for given satellite name
def get_hirs_reader(satellite_name):
    if satellite_name in ['tirosn','noaa06','noaa07','noaa08','noaa09','noaa10','noaa11','noaa12','noaa13','noaa14']:
        print('instrument: HIRS2')
        read_hirs = typhon.datasets.tovs.HIRS2(satname=satellite_name)
    elif satellite_name in ['noaa15','noaa16','noaa17']:
        print('instrument: HIRS3')
        read_hirs = typhon.datasets.tovs.HIRS3(satname=satellite_name)
    elif satellite_name in ['noaa18','noaa19','metopa','metopb','metobc']:
        print('instrument: HIRS4')
        read_hirs = typhon.datasets.tovs.HIRS4(satname=satellite_name)
    else:
        print(f'**** {satellite_name} does not have HIRS on bord ****')
    
    return read_hirs

# calculate the average of counts of given scanline and range 
# of scanpositions as well as the standard deviation
def calc_mean_counts(scanline_counts,timeindex,scan_range=range(10,56)):
    scanline_counts = scanline_counts.isel(time=timeindex)
    mean = np.mean(scanline_counts.sel(scanpos=scan_range).data)
    std = (np.std(scanline_counts.sel(scanpos=scan_range).data)\
        **2/len(scan_range)+0.25/3)**0.5
    return mean, std

# calculate the average of counts of the dsv scanline before and
# after the moon intrusion as well as the standard deviation
def calc_mean_counts_dsv(scanline_counts,timeindex):
    counts_dsv_before_intrusion = scanline_counts.isel(
        time=timeindex - 1).values
    counts_dsv_after_intrusion = scanline_counts.isel(
        time=timeindex + 1).values

    mean = (np.mean(counts_dsv_before_intrusion) + \
        np.mean(counts_dsv_after_intrusion))/2
    std = (((np.std(counts_dsv_before_intrusion)**2 + \
        np.std(counts_dsv_after_intrusion)**2)/48+0.5/3)/2)**0.5
    return mean, std

# calculate the blackbody-target temperature of given timeindex
def calc_bb_temp(data,timeindex,nr_sensors):
    temp_bb = data[{"time":data["scantype"]==1}]\
        ["temperature_internal_warm_calibration_target"].isel(
        time=timeindex).sel(prt_number_iwt=range(0,nr_sensors)).values

    mean = np.mean(temp_bb)
    std = np.std(temp_bb/(nr_sensors+1)**0.5)
    return mean, std

# calculate the blackbody-target flux / radiance of given
# channel and blackbody-target temperature. Correction
# factors are applied to correct the bb temperature.
def calc_bb_flux(ch_wn,corr_factors,bb_temp):
    Temp_corrected = corr_factors[0] + corr_factors[1]*bb_temp
    flux = typhon.physics.planck_wavenumber(ch_wn,Temp_corrected)/2.9979245e8*1E20*0.98
    return flux

# get the satellite longitude and altitude at given timestamp
def get_position(data,timeindex):
    longitude = data[{"time":data["scantype"]==1}]["longitude"].isel(
        time=timeindex).sel(
        scanpos=28).values
    latitude = data[{"time":data["scantype"]==1}]["latitude"].isel(
        time=timeindex).sel(
        scanpos=28).values
    altitude = data[{"time":data["scantype"]==1}]["platform_altitude"].isel(
        time=timeindex).values
    return longitude, latitude, altitude