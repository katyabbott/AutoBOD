### Instructions for producing a netCDF of respiration rates from AutoBOD logfiles

The accompanying scripts are particularly useful if you have a series of different AutoBOD runs that you would like to stitch together for various reasons (i.e., you swapped around BOD bottles to increase the number of respiration rate measurements, or you had to restart the AutoBOD). However, they are still useful for calculating the oxygen concentration from phase angle with salinity correction, as well as generating a netCDF of respiration rates for additional analysis.

What you will need:
- `AutoBOD_processing.py` Python script
- Log files from AutoBOD in a directory
- Metadata CSV that relates the log files into the corresponding cast/Niskin etc.

An example of the metadata CSV:

| Bottle | Bottle_ID | Cast_Niskin | Logfile                   | AutoBOD_start_time_UTC | Depth | Notes |
|--------|-----------|-------------|---------------------------|------------------------|-------|-------|
| 1      | C6N15_2   | C6N15       | Run003_022122_autoBOD.log | 2/21/22 8:08           | 75    |       |
| 2      | C6N15_3   | C6N15       | Run003_022122_autoBOD.log | 2/21/22 8:08           | 75    |       |
| 3      | C8N16_1   | C8N16       | Run003_022122_autoBOD.log | 2/21/22 8:08           | 100   |       |
| 4      | C7N15_1   | C7N15       | Run003_022122_autoBOD.log | 2/21/22 8:08           | 73    |       |
| 5      | C8N16_2   | C8N16       | Run003_022122_autoBOD.log | 2/21/22 8:08           | 100   |       |
| 6      | C7N15_2   | C7N15       | Run003_022122_autoBOD.log | 2/21/22 8:08           | 73    |       |

Here, the "Bottle" field corresponds to the bottle location recorded in the AutoBOD file; it is important to make sure that these match. The "Bottle_ID" represents a unique ID for a BOD bottle. Here, because I took triplicates from each Niskin bottle, I have identified each bottle by the associated CTD cast and Niskin, along with the triplicate number (1, 2 or 3). The "Cast_Niskin" field allows us to group the triplicates together. The "Logfile" field points the program to where it should read the logfile from. The "AutoBOD_start_time_UTC" refers to when an AutoBOD run was started. This is helpful for cross-referencing with the time recorded on the AutoBOD, which usually starts on the 1st of the year. The depth is given in meters, and the notes category is useful for adding notes on a particular bottle to the netCDF.

The Aanderaa operating manual has salinity corrections for measurements that are not made in freshwater. They can be found [here](http://www.argodatamgt.org/content/download/26531/181223/file/Aanderaa_TD218_OperatingManual_OxygenOptode_3830_3835_3930_3975_4130_4175_RevApr2007.pdf) and are applied in the `calc_airsat_o2_conc` function. 

The `autobod_processing` function carries out two steps. First, for a given AutoBOD run, it calculates the oxygen concentration given the temperature and salinity corrections. 

Second, it generates a xarray dataset for each bottle position (up to 12 in total, depending on whether the table is full), which has dimensions corresponding to the number of observations. (For example, we might have a dataset for the bottle located at position 4, with 663 observations of phase, temperature and oxygen concentration). It then pads all of the datasets with NaNs so that they have as many observations as the dataset with the most observations. This allows us to concatenate the datasets into 1 dataset with two dimensions: the number of bottles and the number of observations. 

Then, it starts to build a master dataset of all the observations. For each bottle in the most recent processed run, it checks whether that bottle ID already exists in the master dataset. If it does, then we concatenate the data with the existing data from that bottle ID along the time axis, so that we have a continuous record. If not, then we put that data off to the side for the time being. At the end, we concatenate all data together so that we now have a dataset with the unique bottles (i.e., triplicates for each cast) and number of observations. This allows us to query by individual respiration rates from a triplicate, but also by sets from the same cast.









