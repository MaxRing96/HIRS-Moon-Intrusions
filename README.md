# HIRS-Moon-Intrusions
A collection of Jupyter Notebooks and Python scripts to search HIRS data for moon intrusions.

Note: To be able to run the scripts with their dependencies you should work on thunder.

## 1) hirs_find_moon_intrusions.ipynb
To search given HIRS data for moon intrusions you can use the jupyter notebook hirs_find_moon_intrusions.ipynb.
To do that you specify the satellite name and the channel number as well as the time span for your search in the second code cell of the notebook.
After that you just run the next two code cells and the script will start the search for moon intrusions via a gradient method. 

### Gradient method
To identify possible moon intrusions the notebook loops over all granule files given for the specified time span. For each of these granule files it checks if two conditions are fulfilled. First it checks whether there is a gradient in between Deep-Space-View (DSV) scanlines biger than 50 (photon) counts. If this condition is fullfilled it will look at the counts of the scanline with the minimum average counts. That is the scanline which should include the moon intrusion (intrusion-scanline). The second condition is that the gradients of the DSV counts of the intrusion-scanline should be bigger than zero and below a certain threshold to filter out scanlines with unrealistic varying counts. These are probably instrument errors and can be excluded. Only if these two conditions are met, the notebook will categorize it as a possible moon intrusion.

If there were any possible moon intrusions found in the data, the script will save the paths to the corresponding granule files in a textfile: 'hirs_moon_intrusions_<SATELLITE>_CH<CHANNEL>_<START_DATE>_<END_DATE>.txt'
  
With the last code cell you can plot all the detected moon intrusions. 
  
## 2) hirs_plot_all_channels.py
To further investigate one of the in 1) detected moon intrusions it is useful to see the DSV-counts of the intrusion-scanline in all channels of the satellite. That is what the python script *hirs_plot_all_channels.py* does. It plots a matrix plot of the DSV intrusion-scanline for all channels. The textfile *config_intr_scan_ranges.txt* provides, for every single channel, the range in which is the suspected moon intrusion is visible according to the DSV counts. These ranges are shown as red lines (lower and upper boundary of the range) in the matrix plot and should be manualy adjusted in the textfile for every new moon intrusion. This step is important because also the calculations in 3) are done for these manualy defined ranges!
To run the script type in the command line:
```
python hirs_plot_all_channels.py <PATH-TO-TO-INTRUSION-FILE>
```

## 3) hirs_calc_intrusion_values.py
If the ranges of the moon intrusion are adjusted for every channel in *config_intr_scan_ranges.txt* this script is used to calculate some quantities of the moon intrusion, which are necessary for further investigation. To run the script type in the command line:
```
python hirs_calc_intrusion_values.py <PATH-TO-TO-INTRUSION-FILE>
```
  
