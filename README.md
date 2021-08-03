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
  
  
