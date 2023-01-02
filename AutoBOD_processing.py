import numpy as np
import pandas as pd
import xarray as xr

### Functions

# Function to convert phase angle to oxygen concentration given bottle temperature 
def calc_airsat(phase, temp):
    
    from numpy import tan, pi, sqrt, exp, log
    
    #Defining constants 
    cal0 = 60.21 #B6
    cal100 = 27.42 #B7
    #airpres = 981 #B8
    T0 = 20 #E6
    T100 = 20 #E7
    dF_k = -0.0847 #B12
    f1 = 0.833 #B11
    dksv_k = 0.000416 #B13
    m = 34 #B14
    
    tan_psi0_t100 = tan(((cal0+dF_k*(T100-T0)))*pi/180) #D11
    tan_psi100_t100 = tan(cal100*pi/180) #D13

    A = tan_psi100_t100/tan_psi0_t100*1/m*100**2 #F11
    B = tan_psi100_t100/tan_psi0_t100*100+tan_psi100_t100/tan_psi0_t100*100/m-f1*100/m-100+f1*100 #F12
    C = tan_psi100_t100/tan_psi0_t100-1 #F13

    ksv_t100 = (-B+sqrt(B**2-4*A*C))/(2*A) #H11
    
    air_sat = (-((tan(phase*pi/180))/(tan((cal0+(dF_k*(temp-T0)))*pi/180))*
                        (ksv_t100+(dksv_k*(temp-T100)))+(tan(phase*pi/180))/
                        (tan((cal0+(dF_k*(temp-T0)))*pi/180))*1/m*(ksv_t100+(dksv_k
                        *(temp-T100)))-f1*1/m*(ksv_t100+(dksv_k*(temp-T100)))-(ksv_t100
                        +(dksv_k*(temp-T100)))+f1*(ksv_t100+(dksv_k*(temp-T100))))+(
                        sqrt(((((tan(phase*pi/180))/(tan((cal0+(dF_k*(temp-T0)))*pi/180))*(
                        ksv_t100+(dksv_k*(temp-T100)))+(tan(phase*pi/180))/(tan((cal0+(dF_k*(
                        temp-T0)))*pi/180))*1/m*(ksv_t100+(dksv_k*(temp-T100)))-
                        f1*1/m*(ksv_t100+(dksv_k*(temp-T100)))-(ksv_t100+(dksv_k*
                        (temp-T100)))+f1*(ksv_t100+(dksv_k*(temp-T100))))**2))-
                        4*((tan(phase*pi/180))/(tan((cal0+(dF_k*(temp-T0)))*pi/180))*1/m*
                        ((ksv_t100+(dksv_k*(temp-T100)))**2))*((tan(phase*pi/180))/(
                        tan((cal0+(dF_k*(temp-T0)))*pi/180))-1))))/(2*((tan(phase*pi/180))/
                        (tan((cal0+(dF_k*(temp-T0)))*pi/180))*1/m*((ksv_t100+(dksv_k*(temp-T100)))**2)))
    
    return air_sat

def calc_o2_conc(air_sat, temp, Sal):
    from numpy import exp, log
    #Salinity-related corrections
    # http://www.argodatamgt.org/content/download/26531/181223/file/Aanderaa_TD218_OperatingManual_OxygenOptode_3830_3835_3930_3975_4130_4175_RevApr2007.pdf
    Ts_eq = lambda t: np.log((298.15 - t)/(273.15 + t)) #scaled temperature, depends on temperature of incubation
    B0 = -6.24097e-3
    C0 = -3.1168e-7
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    airpres = 981 #B8


    ## Given the % air sat and temperature of the bottle, we can calculate concentration 
    o2_conc = (((airpres-exp(52.57-6690.9/(273.15+temp)-4.681*log(273.15+temp)))/1013)*
           air_sat/100.*0.2095*(48.998-1.335*temp+0.02755*temp**2-
                                 0.000322*temp**3+0.000001598*temp**4)*32/22.414)
    
    #Salinity correction
    Ts = Ts_eq(temp)
    o2_conc*=exp(Sal*(B0 + B1*Ts+ B2*Ts**2 + B3*Ts**3) + C0*Sal**2)
    
    o2_conc *= 31.25 #umol/L
    
    return o2_conc

def bottle_add_data(ds, log_sheet):
    ds['Depth'] = (('N_bottles',), log_sheet['Depth'])
    notes = log_sheet['Feature'].astype('str').replace({'nan': ' '})
    ds['Feature'] = (('N_bottles',), notes)
    ds['Sample_collection_time'] = (('N_bottles',), log_sheet['Sample_collection_time'])
    ds['SA'] = (('N_bottles',), log_sheet['SA'])
    ds['CT'] = (('N_bottles',), log_sheet['CT'])
    ds['Cast_Niskin'] = (('N_bottles'), log_sheet['Cast_Niskin'])
    UTC_dates = np.repeat(pd.to_datetime(log_sheet['AutoBOD_start_time_UTC']).values, 
          ds.dims['N_obs']).reshape((ds.dims['N_bottles'], 
                    ds.dims['N_obs'])) + ds['Elapsed_time_seconds']
    #ds['UTC_datetime'] = UTC_dates
    return ds

#This function takes in a list of datasets with different number of observations and pads them to the length of
#the most number of obs, then outputs the result as another list 

def standardize_Nobs(ds_list):
    max_nobs = max([ds.dims['N_obs'] for ds in ds_list])
    ds_list_standardized = []
    for ds in ds_list:
        ds = ds.pad(N_obs = (0, max_nobs - ds.dims['N_obs']))
        ds_list_standardized.append(ds)
    return ds_list_standardized, max_nobs

def datetime_from_date_hour(date_hour_df):
    # Input a pandas series with date ('M/D/y') and hour ('H:M:S') date
    # convert into datetime object
    hh_mm_ss = date_hour_df['Hour'].str.split(':', expand=True)
    datetime_df = pd.to_datetime(date_hour_df['Date']) + \
                        pd.to_timedelta(hh_mm_ss[0].astype('int')*60*60 + hh_mm_ss[1].astype('int')*60 + \
                    hh_mm_ss[2].astype('int'), unit='sec')
    return datetime_df

def autobod_processing(log_dir, log_sheet, autoBOD_start_date, error_codes=[1, 5], logging_on=False):
    Sal = 38.4 #adjust per location

    log_sheet_unique = log_sheet[['Bottle_ID', 'Depth', 'AutoBOD_start_time_UTC', 'Feature', 
                    'Cast_Niskin', 'Sample_collection_time','CT', 'SA']].drop_duplicates(subset=['Bottle_ID'])
    
    for i, logfile in enumerate(sorted(list(set(log_sheet['Logfile'].values)))):
        if logging_on:
            print(logfile)

        ### Calculate oxygen concentration and other parameters of interest ###
        
        cols = ['Amplitude', 'Phase', 'Optode_Calculated_O2', 'Error', 'Encoder', 'Bottle', 'Sample', 
            'Date', 'Hour', 'IRDetT', 'IRBotT', 'Steps', 'Light']
        run = pd.read_csv(log_dir + logfile, skiprows=4, skipfooter=3, names = cols,
                    sep='\s+', engine='python')
        run['Phase'] = run['Phase']/100 #Divide phase by 100 to get accurate results

        run['autoBOD_datetime'] = datetime_from_date_hour(run[['Date', 'Hour']])

        # Get relative number of seconds elapsed by subtracting the time logged by the autoBOD from its start time
        run['Elapsed_time_seconds'] = run['autoBOD_datetime'] - autoBOD_start_date

        # Remove all data where error is not in the list of error codes provided (default 1 or 5)
        run = run.where(run.Error.isin(error_codes))
        run = run.dropna()

        # Get a counter of the number of measurements made (i.e., because the autoBOD records 
        # multiple data points when it finds a spot, and we might want to average these measurements). For example, if it takes 10 measurements
        # the measurement number will increment by 1 every 10 observations (when the Sample field resets)
        run['measurement_number'] = ((run['Sample'] - run['Sample'].shift(1)) < 0).cumsum() + 1

        ### Extract data by bottle and put into xarray dataset so we can track auxiliary variables ###
        
        run_log = log_sheet[log_sheet['Logfile'] == logfile] #log_sheet just for this particular run on the AutoBOD

        ds_list = []

        for bottle in run_log['Bottle']:
            df = run[run['Bottle'] == bottle] #select just data from one bottle
            bottle_id = run_log[run_log['Bottle'] == bottle]['Bottle_ID'].values #get bottle ID
            ds = df.to_xarray().drop('index').rename_dims({'index': 'N_obs'}) #convert to xarray with index as the number of measurements
            ds = ds.expand_dims("N_bottles").assign_coords(Bottle_ID=("N_bottles", bottle_id)) # add additional dimension for number of bottles
            # this will make it easier to concatenate down the line
            utc_start_time = run_log[run_log['Bottle'] == bottle]['AutoBOD_start_time_UTC'].values #start time as recorded in the log_file
            ds['UTC_datetime'] = (("N_obs",), pd.to_datetime(utc_start_time[0]) 
                                  + run[run['Bottle'] == bottle]['Elapsed_time_seconds'].values) #adding elapsed time to the recorded start time so that we
            # have the observations in UTC time, so we can concatenate with other runs
            ds_list.append(ds)

        # If this is the first file, then pad to the length of the bottle with the most observations, then concatenate into one dataset    
        if i == 0:
            ds_list, _ = standardize_Nobs(ds_list)
            master_ds = xr.concat(ds_list, dim = 'N_bottles')

        # Else if we already have an existing dataset, we want to concatenate observations from the same bottle together
        elif i > 0:
            updated_ds_list = []
            nonupdated_dataset = master_ds.copy()

            for ds in ds_list:
                bottle_id = ds.Bottle_ID.values
                if bottle_id in master_ds.Bottle_ID.values: # are there previous observations of this bottle_id in the master_ds already?
                    # Remove observations for this bottle from the master dataset, which we have now called nonupdated_dataset
                    nonupdated_dataset = nonupdated_dataset.where(nonupdated_dataset.Bottle_ID != bottle_id, drop=True)
                    updated_ds = xr.concat([master_ds.sel(N_bottles = 
                                master_ds.Bottle_ID == bottle_id), ds], dim='N_obs') # concatenate the new observations for this bottle with the 
                    # existing observations
                    updated_ds_list.append(updated_ds)
                else:
                    updated_ds_list.append(ds)

            updated_ds_list, max_nobs = standardize_Nobs(updated_ds_list) # pad to the length of the bottle with the most observations
            nonupdated_dataset = nonupdated_dataset.pad(N_obs = (0, max_nobs - nonupdated_dataset.dims['N_obs'])) #also padding this dataset so it has 
            # the same number of observations
            master_ds = xr.concat(updated_ds_list + [nonupdated_dataset], dim='N_bottles') #concatenate both datasets

    master_ds = master_ds.sortby('Bottle_ID')
    log_sheet_unique = log_sheet_unique.sort_values(by='Bottle_ID')
    master_ds = bottle_add_data(master_ds, log_sheet_unique) #Add relevant metadata, including depth of the measurement and cast and Niskin
    
    return master_ds





