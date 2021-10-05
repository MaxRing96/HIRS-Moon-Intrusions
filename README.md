# HIRS-Moon-Intrusions
A collection of Jupyter Notebooks and Python scripts to search HIRS data for moon intrusions.

Note: 
- To be able to run the scripts with their dependencies you should work on thunder. 
- The *hirs_src.py* script contains all the functions used in the other scripts.
- The HIRS data can be found on thunder under `/scratch/uni/u237/data/hirs/`.
- For questions concerning the physical and technical theory behind the detection of moon intrusions in HIRS data, please contact Martin Burgdorf (martin.joerg.burgdorf@uni-hamburg.de).
- For questions/problems concerning the code etc. feel free to contact me (Maximilian Ringel - maximilian.ringel@studium.uni-hamburg.de). 

## 1) hirs_find_moon_intrusions.ipynb
To search given HIRS data (the path to the HIRS data is already included in the notebook) for moon intrusions you can use the jupyter notebook *hirs_find_moon_intrusions.ipynb*. To do that you specify the satellite name and the channel number as well as the time span for your search in the second code-cell of the notebook. After that you just run the next two code-cells and the script will start the search for moon intrusions via a gradient method. 

### Gradient method
To identify possible moon intrusions the notebook loops over all HIRS files given for the specified satellite, channel and time span. For each of these granule files it checks if two conditions are fulfilled. First it checks whether there is a gradient in between Deep-Space-View (DSV) scanlines biger than 50 (photon) counts. If this condition is fullfilled it will look at the counts of the scanline with the minimum average counts. That is the scanline which should include the moon intrusion (intrusion-scanline). The second condition is that the gradients of the DSV counts of the intrusion-scanline should be bigger than zero and below a certain threshold to filter out scanlines with unrealistic varying counts. These are probably instrument errors and can be excluded. Only if these two conditions are met, the notebook will categorize it as a possible moon intrusion.

If there were any possible moon intrusions found in the data, the script will save the paths to the corresponding HIRS files in a textfile: 'hirs_moon_intrusions_<SATELLITE>_CH<CHANNEL>_<START_DATE>_<END_DATE>.txt'
  
With the last code-cell you can plot all the detected moon intrusions as plots like this example:

<img width="1152" alt="Screenshot 2021-08-03 at 13 02 17" src="https://user-images.githubusercontent.com/62293752/128007264-8a6fe7b7-cdee-46c7-8720-5b6ff51267fd.png">
  
All the plots of the moon intrusions are saved under the directory *hirs_moon_intrusion_suspects_plots/* which will be created in your local directory if it doesn't exist already.

To further investigate one of the detected intrusions, you just need the corresponding path to the HIRS file (PATH-TO-TO-INTRUSION-FILE) to continue with step 2)- 4).
  
## 2) hirs_plot_all_channels.py
It is useful to see the DSV-counts of the intrusion-scanline in all channels of the satellite. That is what the python script *hirs_plot_all_channels.py* does. It plots a matrix plot of the DSV intrusion-scanline for all channels (see example below). The textfile *config_intr_scan_ranges.txt* provides, for every single channel, the range in which is the suspected moon intrusion is visible according to the DSV counts. These ranges are shown as red lines (lower and upper boundary of the range) in the matrix plot and should be manualy adjusted in the textfile for every new moon intrusion. The complete scan range reaches from scanposition 10 to scanposition 56. The upper and lower boundary of the intrusion range can be selected in between this range. This step is important because also the calculations in 3) are done for these manualy defined ranges! The output plot is saved under the name *hirs_moon_intrusion_<SATELLITE>_<DATE_OF_INTRUSION>_<SCAN_TIME>_allchannels_dsv.pdf*
To run the script type in the command line:
```
python hirs_plot_all_channels.py <PATH-TO-TO-INTRUSION-FILE>
```

An example plot is shown here:
  
<img width="827" alt="Screenshot 2021-08-03 at 13 18 02" src="https://user-images.githubusercontent.com/62293752/128007461-29949e16-a0da-4780-9fe3-a3e4aa9c55e2.png">

## 3) HORIZON Web-page
Next step is to go to Horizon webpage https://ssd.jpl.nasa.gov/horizons.cgi and get the ANGULAR DIAMETER and the PHASE ANGLE of the moon.
You will need the following input variables:
```
Ephemeris Type: OBSERVER
Target Body: Moon[Luna][301]
Observer Location: Topocentric (lon, lat, alt)
Time span: Date of intrusion, Date_of_intrusion+1min, Step=1minute')
Table Settings: QUANTITIES= 13 (Target angular diamenter), 23 (Sun-Observer-Target), 24 (Sun-Target-Observer ~PHASE angle)
Display/Output: default
```
<img width="1408" alt="Bildschirmfoto 2021-09-28 um 14 53 23" src="https://user-images.githubusercontent.com/90314017/136040330-98948e44-3eb3-4796-8246-fd7a12fa0b69.png">
<img width="947" alt="Bildschirmfoto 2021-09-28 um 14 53 10" src="https://user-images.githubusercontent.com/90314017/136040352-9bf1aac0-b1eb-4a50-842f-ba5f16140d2a.png">

## 4) hirs_calc_intrusion_values.py
If the ranges of the moon intrusion are adjusted for every channel in *config_intr_scan_ranges.txt*, this script can be used to calculate the quantities of the moon intrusion. For the calculation of the radiance of the blackbody-target of the satellite, additional information is needed - The central wavelengths and correction coefficients for every channel and every satellite. These are provided in the folder *meta_files*. If there is no such meta file provided for the choosen satellite, the script will raise a warning but will still calculate all the other values. Only the blackbody radiance will be left out then. 
For the calculation of the radiance and the brightness temperature of the moon, the FOV and the angular diameter is needed. Therefore the script has multiple inputs. Besides the path to the intrusion file, the angular diameter in arcsec is an input as well as the phase-angle of the moon in degree with appended information, whether the moon in Trailing or Leading the Sun (T/L).
The calculated values are saved to a csv file with the name *hirs_moon_intrusion_<SATELLITE>_<DATE_OF_INTRUSION>_<INTRUSION_TIMESTAMP>_calculations.csv* (see example below). To run the script type in the command line:
```
python hirs_calc_intrusion_values.py <PATH-TO-TO-INTRUSION-FILE> <ANG_DIAM> <PHASE-ANGLE>
Example: python hirs_calc_intrusion_values.py /scratch/uni/u237/data/hirs/noaa15_hirs_2020/02/09/NSS.HIRX.NK.D20040.S0111.E0256.B1308081.GC.gz 1961.271 5.4418T
```
The csv table containing the calculated values for the moon intrusion looks like this:
<img width="1873" alt="Bildschirmfoto 2021-10-05 um 16 16 25" src="https://user-images.githubusercontent.com/90314017/136042808-543a3a5f-077d-43f2-ab65-ec88b0326d0c.png">

   
