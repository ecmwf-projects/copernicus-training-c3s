#!/usr/bin/env python
# coding: utf-8

# ![logo](./img/LogoLine_horizon_C3S.png)

# <br>

# # Seasonal Forecast Anomalies

# ### About

# This notebook provides a practical introduction to calculating seasonal forecast anomalies with data from the Copernicus Climate Change Service (C3S). C3S seasonal forecast products are based on data from several state-of-the-art seasonal prediction systems. In this tutorial we shall focus on the [ECMWF SEAS5 model](https://confluence.ecmwf.int/display/CKB/Description+of+SEAS5+C3S+contribution), which is one of the forecasting systems available through C3S.
# 
# The tutorial will demonstrate how to access real-time forecast data of total precipitation, with a forecast start date in May 2021 and 6 monthly lead times (up to October 2021). Hindcast data for the same start date and lead-time months in the reference period 1993 to 2016 will also be downloaded. The tutorial will then show how to extract a subset area over South Asia for both the forecast and hindcast data. The climate mean for the reference period will be computed and this reference mean will be subtracted from the real-time forecast data to derive monthly anomalies. These will be visualised as both spatial maps and time series. Finally, 3-monthly anomalies will be calculated and visualised in an interactive plot, as a demonstration of how to reproduce similar [charts available through C3S](https://climate.copernicus.eu/charts/c3s_seasonal/).

# 
# The notebook has the following outline:
# * 1 - Download data from the CDS
# * 2 - Hindcast data processing: calculate the reference climate mean
# * 3 - Real-time forecasts: calculate seasonal forecast anomalies
# * 4 - Visualize seasonal forecast monthly anomalies for a geographical subregion
#   * 4.1 - Spatial maps
#   * 4.2 - Time series of regional averages
# * 5 - Reproduce C3S graphical products: compute 3-month anomalies

# <br>

# <div class="alert alert-block alert-success">
# <b>NOTE</b>: <br>
#     <a href="https://cds.climate.copernicus.eu/cdsapp#!/dataset/seasonal-postprocessed-single-levels">Precomputed anomalies are also available through the CDS</a>. Note these may be slightly different due to minor differences in the way they are computed (e.g. months of constant length, 30 days)  and also due to GRIB packing discretisation. <a href="https://confluence.ecmwf.int/display/UDOC/Why+are+there+sometimes+small+negative+precipitation+accumulations+-+ecCodes+GRIB+FAQ">See here for more detials.</a></div>

# Please see here the full documentation of the [C3S Seasonal Forecast Datasets](https://confluence.ecmwf.int/display/CKB/C3S+Seasonal+Forecasts%3A+datasets+documentation). This notebook introduces you to the [seasonal forecast monthly statistics](https://cds.climate.copernicus.eu/cdsapp#!/dataset/seasonal-monthly-single-levels?tab=overview) datasets on single levels (as opposed to multiple levels in the atmosphere).

# ### How to access the notebook
# 
# This tutorial is in the form of a [Jupyter notebook](https://jupyter.org/). You will not need to install any software for the training as there are a number of free cloud-based services to create, edit, run and export Jupyter notebooks such as this. Here are some suggestions (simply click on one of the links below to run the notebook):

# |Binder|Kaggle|Colab|NBViewer|
# |:-:|:-:|:-:|:-:|
# |[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ecmwf-projects/copernicus-training/HEAD?urlpath=lab/tree/C3S_seasonal-forecasts-anomalies.ipynb)|[![Kaggle](https://kaggle.com/static/images/open-in-kaggle.svg)](https://kaggle.com/kernels/welcome?src=https://github.com/ecmwf-projects/copernicus-training/blob/master/C3S_seasonal-forecasts-anomalies.ipynb)|[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ecmwf-projects/copernicus-training/blob/master/C3S_seasonal-forecasts-anomalies.ipynb)|[![NBViewer](https://raw.githubusercontent.com/ecmwf-projects/copernicus-training/master/img/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/ecmwf-projects/copernicus-training/blob/master/C3S_seasonal-forecasts-anomalies.ipynb)|
# |Binder may take some time to load, so please be patient!|You will need to login/register, and switch on the internet via *settings*|You will need to run the command `!pip install cartopy` before importing the libraries|This will not run the notebook, only render it|

# If you would like to run this notebook in your own environment, we suggest you install [Anaconda](https://docs.anaconda.com/anaconda/install/), which contains most of the libraries you will need. You will also need to install [Xarray](http://xarray.pydata.org/en/stable/) for working with multidimensional data in netcdf files, and the CDS API (`pip install cdsapi`) for downloading data programatically from the CDS.

# <hr>

# ### Install packages

# In[2]:


# Install CDS API for downloading data from the CDS
get_ipython().system('pip install cdsapi')


# In[1]:


# Install cfgrib to enable us to read GRIB format files
get_ipython().system('conda install -c conda-forge cfgrib -y')


# ### Load packages

# In[45]:


# Miscellaneous operating system interfaces
import os

# CDS API
import cdsapi

# To map GRIB files to NetCDF Common Data Model
import cfgrib

# Libraries for working with multi-dimensional arrays
import numpy as np
import xarray as xr
import pandas as pd

# Libraries for plotting and geospatial data visualisation
import matplotlib.path as mpath
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature

# To work with data labels in dictionary format
from collections import OrderedDict

# Date and time related libraries
from dateutil.relativedelta import relativedelta
from calendar import monthrange
import datetime

# Interactive HTML widgets
import ipywidgets as widgets

# Disable warnings for data download via API
import urllib3 
urllib3.disable_warnings()


# <hr>

# ## 1. Request data from the CDS programmatically with the CDS API

# The first step is to request data from the Climate Data Store programmatically with the help of the CDS API. Let us make use of the option to manually set the CDS API credentials. First, you have to define two variables: `URL` and `KEY` which build together your CDS API key. Below, you have to replace the `#########` with your personal CDS key. Please find [here](https://cds.climate.copernicus.eu/api-how-to) your personal CDS key.

# In[ ]:


URL = 'https://cds.climate.copernicus.eu/api/v2'
KEY = '########################################'


# Here we specify a data directory in which we will download our data and all output files that we will generate:

# In[62]:


DATADIR = '../input/seasonal'


# The next step is then to request the seasonal forecast monthly statistics data on single levels with the help of the CDS API. Below, we download two separate files of total precipitation for six monthly lead times (start date in May):
# * **Retrospective forecasts (Hindcasts) for 1993 to 2016**
# * **Forecasts for 2021**
# 
# Seasonal forecast data are disseminated in the GRIB data format. Let us store the data in the main working directory with the names: 
# * `ecmwf_seas5_1993-2016_05_hindcast_monthly_tp.grib` and 
# * `ecmwf_seas5_2021_05_forecast_monthly_tp.grib`.
# 
# Running the code block below will download the data from the CDS as specified by the following API keywords:
# 
# > **Format**: `Grib` <br>
# > **Originating centre**: `ECMWF` <br>
# > **System**: `5` *this refers to SEAS5* <br>
# > **Variable**: `Total precipitation` <br>
# > **Product type**: `Monthly mean` *all ensemble members will be retrieved* <br>
# > **Year**: `1993 to 2016` *for the hindcast* `2021` *for the forecast* <br>
# > **Month**: `05` *May* <br>
# > **Leadtime month**: `1 to 6` *May to October*

# If you have not already done so, you will need to accept the **terms & conditions** of the data before you can download it. These can be viewed and accepted in the [CDS download page](https://cds.climate.copernicus.eu/cdsapp#!/dataset/seasonal-monthly-single-levels?tab=form) by scrolling to the end of the download form.

# <div class="alert alert-block alert-success">
# <b>NOTE</b>: <br>
#     The API request below can be generated automatically from the <a href="https://cds.climate.copernicus.eu/cdsapp#!/dataset/seasonal-monthly-single-levels?tab=form">CDS download page</a>. At the end of the download form there is a <code>Show API request</code> icon, which allows a copy-paste of the code below.</div>

# In[ ]:


c = cdsapi.Client(url=URL, key=KEY)

# Hindcast data request
c.retrieve(
    'seasonal-monthly-single-levels',
    {
        'format': 'grib',
        'originating_centre': 'ecmwf',
        'system': '5',
        'variable': 'total_precipitation',
        'product_type': 'monthly_mean',
        'year': [
            '1993', '1994', '1995',
            '1996', '1997', '1998',
            '1999', '2000', '2001',
            '2002', '2003', '2004',
            '2005', '2006', '2007',
            '2008', '2009', '2010',
            '2011', '2012', '2013',
            '2014', '2015', '2016',
        ],
        'month': '05',
        'leadtime_month': [
            '1', '2', '3',
            '4', '5', '6',
        ],
    },
    f'{DATADIR}/ecmwf_seas5_1993-2016_05_hindcast_monthly_tp.grib')

# Forecast data request
c.retrieve(
    'seasonal-monthly-single-levels',
    {
        'format': 'grib',
        'originating_centre': 'ecmwf',
        'system': '5',
        'variable': 'total_precipitation',
        'product_type': 'monthly_mean',
        'year': '2021',
        'month': '05',
        'leadtime_month': [
            '1', '2', '3',
            '4', '5', '6',
        ],
    },
    f'{DATADIR}/ecmwf_seas5_2021_05_forecast_monthly_tp.grib')


# <br>

# ## 2. Calculate seasonal hindcast climate mean

# Seasonal forecasts are affected by systematic errors (biases) which are dependent on the leadtime, the time of year, the variable and the location. Hindcast data can help us to understand and account for these biases. We will calculate the hindcast climate mean for each lead time month, averaged over the years 1993 to 2016. In the next section we will then calculate the anomalies, i.e. the deviation of the 2021 forecast for each lead time month with respect to the hindcast mean. Anomalies are thus calculated "in the model space" where the model climate calculated from the hindcasts is taken as the reference.
# 
# #### Read the downloaded data
# 
# We will use the Python library [xarray](http://xarray.pydata.org/en/stable/) and its function `open_dataset` to read the GRIB file of the hindcast data, specifying the keyword argument `engine` and `'cfgrib'`. [cfgrib](https://github.com/ecmwf/cfgrib) is a Python interface to map GRIB files to the NetCDF Common Data Model using [ecCodes](https://github.com/ecmwf/eccodes).
# 
# The result is a `xarray.Dataset` object with five dimensions:
# 
# > **Number**: Ensemble members (25) <br>
# > **Time**: Forecast start date for each year (1st of May) <br>
# > **Step**: Lead time (nanoseconds in each leadtime month) <br>
# > **Latitude**: Latitudes in 1 deg resolution<br>
# > **Longitude**: Longitudes in 1 deg resolution and in 0-360 grid<br>

# In[63]:


ds = xr.open_dataset(f'{DATADIR}/ecmwf_seas5_1993-2016_05_hindcast_monthly_tp.grib', engine='cfgrib')
ds


# #### Change representation of forecast lead time
# 
# Notice that xarray attempts to extract the time dimensions from the metadata in the GRIB file. By default, it tries to arrange the shape of the data array using “time” (the start time) and “step” (the lead-time in nanoseconds). This comes with two potential issues:
# - Due to some limitations in the GRIB edition 1 encoding of long forecasts, monthly aggregations such as these do not perfectly encode the start and end of the aggregation period (calendar months), and instead they use “step”, which points to the end of the aggregating interval. This might be misleading.  
# - Due to the different lengths of February in leap and non/leap years, we might have different values of “step” for the same “leadtime_month” index when we put together different start years, as we do here.
# 
# A more useful representation of the data is to replace the dimension “step” with “forecastMonth”, which is part of the GRIB metadata and provides the same integer indices used in the CDS API syntax for “leadtime_month”, with the value 1 corresponding to the first complete calendar month after the forecast start date. For instance for a forecast with start date on the 1st May, `forecastMonth=1` is May, and for a forecast with start date on 17th April, `forecastMonth=1` would also be May. This is more coherent and avoids the ambiguity described above.
# 
# We will create this custom data structure through use of the  keyword argument `backend_kwargs` with the `time_dims` option, specifying the two time dimensions: `forecastMonth` and `time` (or `indexing_time`, see note below). The new dataset will then have `forecastMonth` instead of `step`, with values from 1 to 6 (for each lead time month). In this conversion we will lose some information about `valid_time`, which we will need to calculate later.

# <div class="alert alert-block alert-success">
# <b>NOTE</b>: <br>
#     The second of the time dimensions is valid for systems with burst start dates (such as in our example), but for lagged systems, <i>time</i> should be replaced with <i>indexing_time</i>. Please see <a href="https://confluence.ecmwf.int/display/CKB/Seasonal+forecasts+and+the+Copernicus+Climate+Change+Service#heading-Burstvslaggedmode">here</a> for more details on the difference between burst and lagged systems.</div>

# In[64]:


ds_hindcast = xr.open_dataset(f'{DATADIR}/ecmwf_seas5_1993-2016_05_hindcast_monthly_tp.grib', engine='cfgrib', backend_kwargs=dict(time_dims=('forecastMonth', 'time')))
ds_hindcast


# #### Extract data array from dataset
# 
# The `xarray.Dataset` object into which we have read the data from the downloaded GRIB files may contain arrays of multiple variables (even if we have only one: total precipitation). Another xarray data structure, `xarray.DataArray`, facilitates operations on single variables. We will use this to further process our data. You can select the relevant DataArray from a Dataset by specifying the name of the variable (in our case `tprate`) in square brackets `[]`. 

# In[65]:


tprate_hindcast = ds_hindcast['tprate']
tprate_hindcast


# A DataArray provides additional attributes of the variable. For example, you see that the precipitation is expressed as a rate with the unit m/s. Some variables (like precipitation, radiation or heat fluxes) are encoded in GRIB as accumulations from the beginning of the forecast, and their monthly aggregations are therefore expressed as rates. Later in this tutorial we shall convert this to total accumulation.

# #### Average over ensemble members and years to create hindcast climatology
# 
# We can now create the hindcast climatology by averaging over the 25 ensemble members and the 24 years. We do this for each forecast lead time and for each geographical grid point. We use the function `mean()` to average over one or more dimensions, which in this case are `number` (25 ensemble members) and `time` (years from 1993 to 2016). The result is an `xarray.DataArray` with three dimensions: `forecastMonth`, `latitude` and `longitude`.

# In[66]:


tprate_hindcast_mean = tprate_hindcast.mean(['number', 'time'])
tprate_hindcast_mean


# <br>

# ## 3. Load monthly seasonal forecast data for 2021 and calculate seasonal forecast anomalies

# The next step is to load the real-time seasonal forecast data for 6 lead time months in 2021, beginning in May. We will then subtract the hindcast climatology from this to derive the anomalies. 
# 
# #### Load seasonal forecast data for 2021 and change representation of forecast lead time
# 
# Before we can subtract the hindcast anomaly from the 2021 forecast, we need first to ensure the two datasets have the same structure. We thus need to apply the same processing on the forecast data as we did with the hindcast to use `forecastMonth` as a coordinate of the array instead of `step`.

# In[67]:


seas5_forecast = xr.open_dataset(f'{DATADIR}/ecmwf_seas5_2021_05_forecast_monthly_tp.grib', engine='cfgrib', 
                                 backend_kwargs=dict(time_dims=('forecastMonth', 'time')))
seas5_forecast


# Once the `xarray.Dataset` is loaded, you see that it has 4 dimensions:
# > **Number**: 51 ensemble members <br>
# > **ForecastMonth**: 6 leadtime months <br>
# > **Latitude**: latitude values <br>
# > **Longitude**: longitude values <br>
# 
# Compared to the hindcast data, we have only a single start date (May 2021) in the `time` coordinate. Another difference is that seasonal forecast real-time data have 51 ensemble members, compared to hindcasts, which only have 25 ensemble members.
# 
# #### Calculate 2021 anomalies
# 
# Now, we can compute the seasonal forecast anomalies for the 6 lead time months beginning with May 2021. We can compute the anomalies by subtracting the long term average (`tprate_hindcast_mean`) from the seasonal forecast real-time data for May 2021 (`seas5_forecast`).
# 
# The resulting `xarray.DataArray` has the anomaly information for each of the lead time months from May 2021 to October 2021 relative to the reference period May to October for the years 1993 to 2016.

# In[68]:


seas5_anomalies_202105 = seas5_forecast['tprate'] - tprate_hindcast_mean
seas5_anomalies_202105


# #### Convert forecast lead time month into dates
# 
# A closer look at the coordinate information of `forecastMonth` reveals that it only gives integer values from 1 to 6, but not the actual time for which the forecast is valid. We will need this latter information when we convert from precipitation rate to accumulation, as for this we will need to know the number of days per month. Additionally, this information will be useful for later labelling.
# 
# To calculate the time for which the forecast is valid, let us assign a new coordinate with the name `valid_time` based on the dimension `forecastMonth`. This coordinate information shall provide us the valid timestamp for each of the forecast months. We can use the value of the time dimension as the start date. The function `relativedelta()` allows us to add a specific number of months to the start date.
# 
# In a second step, with the xarray function `assign_coords()`, we assign this newly created `DateTimeIndex` object as a new coordinate.

# In[69]:


valid_time = [pd.to_datetime(seas5_anomalies_202105.time.values) + relativedelta(months=fcmonth-1) for fcmonth in seas5_anomalies_202105.forecastMonth]
seas5_anomalies_202105 = seas5_anomalies_202105.assign_coords(valid_time=('forecastMonth',valid_time))
seas5_anomalies_202105


# #### Convert from precipitation rates to accumulation
# 
# Before visualizing the seasonal forecast anomalies, let us convert the precipitation from rates in m/s to total accumulated precipitation per month in mm. The conversion has to be done per month, as it is dependent on the number of days of a specific month. For this, we create a new dimension called `numdays`. The function `monthrange()` from the calendar library returns the number of days in a month, based on a given year and month. With the help of the new coordinate information `valid_time`, we can thus create a list of the number of days. Again, with the function `assign_coords()`, we assign this list as a new coordinate in our data array.

# In[70]:


numdays = [monthrange(dd.year,dd.month)[1] for dd in valid_time]
seas5_anomalies_202105 = seas5_anomalies_202105.assign_coords(numdays=('forecastMonth',numdays))
seas5_anomalies_202105


# We can now use the newly created coordinate information to convert for each month the precipitation accumulations to total precipitation. For this, we multiply the precipitation values with the number of days and then with 24 x 60 x 60 (number of seconds in a day) to retrieve precipitation values in m. As a last step, we convert the values in m to mm by multiplying them by 1000.

# In[71]:


seas5_anomalies_202105_tp = seas5_anomalies_202105 * seas5_anomalies_202105.numdays * 24 * 60 * 60 * 1000
seas5_anomalies_202105_tp


# As a last step before visualizing the data, we want to add the attributes units and long_name and update them, as they have changed with the previous workflow steps. You can simply specify the name of an attribute (e.g. units ) and assign it a new value.

# In[72]:


seas5_anomalies_202105_tp.attrs['units'] = 'mm'
seas5_anomalies_202105_tp.attrs['long_name'] = 'Total precipitation anomaly' 
seas5_anomalies_202105_tp


# ## 4. Visualize seasonal forecast anomalies for a geographical subregion

# To visualise the seasonal forecast anomalies, we could simply look at the ensemble mean to summarise the information given by all individual ensemble members. However, this would mean we lose the richness of information the whole ensemble provides. To highlight the relevance of considering not only the ensemble mean but the complete distribution, we will plot the anomalies for all individual members. We will do this in two ways:
# 
# > 1. **Spatial map** of all individual members for a given lead time. <br>
# > 2. **Time series plot** of all members averaged over a given region. <br>
# 
# We will explore these two methods in the respective subsections below, but we will do so for a geographical subset for South Asia with a bounding box of `[N:30, W:70, S:5, E:90]`. The first step is therefore to subset the data over this area.
# 
# #### Define a geographical subset
# 
# To create a subset, we will define a function called `ds_latlon_subet()` which generates a geographical subset of an xarray data array. The latitude and longitude values that are outside the defined area are dropped.

# In[73]:


def ds_latlon_subset(ds,area,latname='latitude',lonname='longitude'):
    mask = (ds[latname]<=area[0]) & (ds[latname]>=area[2]) & (ds[lonname]<=area[3]%360) & (ds[lonname]>=area[1]%360)
    dsout = ds.where(mask,drop=True)
    return dsout


# Now, we apply the function `ds_latlon_subset()` to the data array of the seasonal forecast anomalies.

# In[74]:


sub = (30, 70, 5, 90) # North/West/South/East
seas5_SAsia = ds_latlon_subset(seas5_anomalies_202105_tp, sub)
seas5_SAsia


# ### 4.1 Spatial map visualization

# Now we can create a spatial map visualisation which shows for a given lead time the total precipitation anomaly of the 51 ensemble members.

# In[75]:


# Select a leadtime to visualise
lead_time = 2


# The plotting code below can be split into the following sections:
# > 1. **Define the overall figure**: Including the size and spacing between subplots. <br>
# > 2. **Create the plots**: Loop over each subplot and set a customized title, add coastlines and set geographical extent. <br>
# > 3. **Create a colour bar**: Add a colour bar with a label. <br>
# > 4. **Save the figure**: Save the figure as a png file. <br>

# In[104]:


# Define figure and spacing between subplots
fig = plt.figure(figsize=(16, 12))
plt.subplots_adjust(hspace=0.15, wspace = 0.05)

new_date_format = seas5_SAsia['valid_time'].dt.strftime('%b, %Y')

# Define overall title
plt.suptitle('C3S ECMWF SEAS5 total precipitation anomaly'
             + os.linesep + 
             f'Start data: {str(new_date_format[0].data)}' 
             + f' - Valid date: {str(new_date_format[lead_time-1].data)}'
             , fontsize=18)

# Define each subplot looping through the ensemble members
for n in np.arange(51):
    # Add a new subplot iteratively
    ax = plt.subplot(6, 10, n + 1, projection=ccrs.PlateCarree())
    # Plot data
    im = ax.pcolormesh(seas5_SAsia.longitude.values, seas5_SAsia.latitude.values, 
                       seas5_SAsia[n,lead_time-1,:,:], cmap='bwr_r', vmin=-350, vmax=350)
    ax.set_title(f'Member {n+1}') # Set subplot title
    ax.coastlines(color='black') # Add coastlines
    # Set the extent (x0, x1, y0, y1) of the map in the given coordinate system.
    ax.set_extent([sub[1],sub[3],sub[2],sub[0]], crs=ccrs.PlateCarree())

# Create a colour bar at the bottom of the fugure
fig.subplots_adjust(bottom=0.0)
# Add axis to make space for colour bar (left, bottom, width, height)
cbar_ax = fig.add_axes([0.3, 0.05, 0.5, 0.02])
fig.colorbar(im, cax=cbar_ax, orientation='horizontal', label='Total precipitation anomaly (mm)')

# Save the figure
fig.savefig('./TotalPrecAnomalyForecastSAsia.png')


# <br>

# ### 4.2 Plot of total precipitation anomalies for each seasonal forecast month

# In this step we will summarise the total precipitation behaviour over the whole region for each lead time month. We will do this by averaging in the spatial (latitude and longitude) dimensions. 
# 
# To put the anomalies in context they will be compared to the reference climate computed in this subregion from the hindcast data.

# #### Average over region

# The first step is to create the weighted average of the subregion defined above. For this, we first estimate the cell area with the cosine of the latitude. These weights are then applied when the data array is averaged over the two dimensions `latitude` and `longitude`. You can use the xarray function `weighted()` together with the function `mean()` to create a weighted average of the geographical region. The result is a data array with two dimensions: `number` and `forecastMonth`.

# In[77]:


weights =np.cos(np.deg2rad(seas5_SAsia.latitude))

anoms_SAsia = seas5_SAsia.weighted(weights).mean(['latitude','longitude'])
anoms_SAsia


# #### Conversion of data into two dimensional pivot table

# To facilitate further processing we will convert the data array into a pandas.Dataframe object with the function `to_dataframe()`. Pandas dataframes are optimised for the processing and visualisation of two dimensional data. We may want to drop coordinates that are not needed. We can do this with the function `drop_vars()`.

# In[97]:


anoms_SAsia_df = anoms_SAsia.drop_vars(['time','surface','numdays']).to_dataframe(name='anomaly')
anoms_SAsia_df


# We will now reorganise the data structure to facilitate plotting. We will do this through the use of the following Pandas functions:
# 
# > **reset_index().drop()**: We use `reset_index()` to convert the ensemble members into columns, and we use `drop()` to remove the *forecastMonth* column, as we have the same information in the *valid_time* dimension <br>
# > **set_index()**: allows us to define which column(s) shall be used as data frame indices <br>
# > **unstack()**: converts the columns into a pivot table, thus re-arranging the data into an easily manageable two dimensional table.<br>

# In[98]:


anoms_SAsia_df = anoms_SAsia_df.reset_index().drop('forecastMonth',axis=1).set_index(['valid_time','number']).unstack()
anoms_SAsia_df


# The data frame above provides us for each ensemble member and forecast month the total precipitation anomaly.
# 
# As a last step, to improve the labels of our data plot later on, we will convert the `valid_time` column values from *YYYY-MM-DD* format into *Month, Year* format.

# In[102]:


anoms_SAsia_m_yr = anoms_SAsia_df.reset_index()
anoms_SAsia_m_yr['valid_time'] = anoms_SAsia_m_yr['valid_time'].dt.strftime('%b, %Y')
anoms_SAsia_m_yr


# #### Repeat steps for hindcast data (to compare forecast and hindcast anomalies)
# 
# We would like now to repeat the steps for the hindcast data in order to show the forecast in the context of the reference climate. I.e. we will compare the seasonal forecast and hindcast anomalies for each lead-time month. In the final plot we will show the anomaly of each forecast ensemble member, while the model climate defined by the hindcasts will be represented in the form of tercile categories. This will put each ensemble member forecast anomaly into context and will give us an indication of its intensity.
# 
# We therefore revisit the hindcast data (`tprate_hindcast`), which we loaded in [section 2](#hindcast_climate_mean) above. We will apply the same steps to create a geographical subset and aggregate spatially.

# In[24]:


tprate_hindcast_SAsia = ds_latlon_subset(tprate_hindcast, sub)


# To create the weighted average of the two dimensions `latitude` and `longitude`, we can use the same `weights` as above.

# In[25]:


tprate_hindcast_SAsia = tprate_hindcast_SAsia.weighted(weights).mean(['latitude','longitude'])


# The resulting array has three remaining dimensions, `number`, `forecastMonth` and `time`.

# The next step is to compute the hindcast anomalies using the hindcast climatology as reference. As before, we create the hindcast climatology by computing the mean over the two dimensions `number` (ensemble member) and `time` (years) resulting in the climatology for each lead-time month. We then subtract this from the hindcast data to derive the anomalies.

# In[26]:


anom_hindcast = tprate_hindcast_SAsia - tprate_hindcast_SAsia.mean(['number','time'])


# We now repeat the same steps to convert precipitation accumulations in m/s to total precipitation in m.

# In[27]:


valid_time = [ pd.to_datetime(anom_hindcast.time.values[0]) + relativedelta(months=fcmonth-1) for fcmonth in anom_hindcast.forecastMonth]
numdays = [monthrange(dd.year,dd.month)[1] for dd in valid_time]
anom_hindcast = anom_hindcast.assign_coords(numdays=('forecastMonth',numdays))


# In[28]:


anom_hindcast_tp = anom_hindcast * anom_hindcast.numdays * 24 * 60 * 60 * 1000


# #### Calculation of hindcast anomaly terciles
# 
# Now, for each of the lead-time months, we can compute the minimum and maximum as well as the 1/3 and 2/3 percentiles of the hindcast anomaly ensemble.

# In[29]:


P0 = anom_hindcast_tp.min(['number','time'])
P33 = anom_hindcast_tp.quantile(1/3.,['number','time'])
P66 = anom_hindcast_tp.quantile(2/3.,['number','time'])
P100 = anom_hindcast_tp.max(['number','time'])


# A last step before we visualise this data is to 

# The last step is to visualize the hindcast climatology boundaries as filled area and the ensemble seasonal forecast anomalies for each forecast month. The visualisation code below has three main sections:
# > 1. **Initiate the figure**: Initiate a matplotlib figure <br>
# > 2. **Plot the seasonal forecast anomalies + climatology boundaries**: Plot the ensemble seasonal forecast anomalies as black circles and the hindcast climatology as filled areas <br>
# > 3. **Customise plot settings**: Add additonal features, such as title, x- and y-axis labels etc. <br>

# In[103]:


# Initiate the figure
fig=plt.figure(figsize=(20, 10))
ax=fig.gca()

# Plot the seasonal forecast anomalies + climatology boundaries as filled areas
ax.plot(anoms_SAsia_m_yr.valid_time, anoms_SAsia_m_yr.anomaly, marker='o', linestyle='', color='black', alpha=0.5, label='Forecast')
ax.fill_between(anoms_SAsia_m_yr.index,P66,P100,color='#5ab4ac',label='Above normal')
ax.fill_between(anoms_SAsia_m_yr.index,P33,P66,color='lightgray',label='Near average')
ax.fill_between(anoms_SAsia_m_yr.index,P0,P33,color='#d8b365',label='Below normal')

# Customize plot settings
plt.title('C3S ECMWF SEAS5 total precipitation anomaly (South Asia)' 
          + os.linesep + 
          'Start date: May 2021 - Reference period: 1993-2016' 
          + os.linesep, fontsize=14)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

ax.set_ylabel('Total precipitation anomaly (mm)' + os.linesep, fontsize=12)
ax.set_xlabel(os.linesep + 'Valid date', fontsize=12)

# Save the figure
fig.savefig('./TotalPrecForecastHindcastAnomaliesSAsia.png')


# <br>

# ## 5. Calculate seasonal forecast 3-month anomalies

# In this last section of the tutorial we will show how to create plots similar to those on the [Copernicus Climate Change Service (C3S) website](https://climate.copernicus.eu/charts/c3s_seasonal/). This process includes calculating 3-month aggregations and ensemble means. This differs from the previous sections where we have focussed on the processing of monthly data and individual members.
# 
# #### Compute 3-month rolling averages
# 
# The first step is to compute 3-month rolling averaged for the seasonal forecasts and hindcast data. We will do this through a combination of the xarray functions `rolling()` and `mean()`.

# In[33]:


seas5_forecast_3m = seas5_forecast.rolling(forecastMonth=3).mean()
ds_hindcast_3m = ds_hindcast.rolling(forecastMonth=3).mean()


# #### Calculate anomalies
# 
# We now repeat the process of calculating anomalies with respect to the climatology, this time for the 3-monthly data.

# In[34]:


ds_hindcast_3m_hindcast_mean = ds_hindcast_3m.mean(['number','time'])
seas5_anomalies_3m_202105 = seas5_forecast_3m.tprate - ds_hindcast_3m_hindcast_mean.tprate


# #### Ensemble mean anomaly
# 
# We want to compute the average of the 3-month seasonal forecast anomaly of the ensemble members. For this, we have to average over the dimension `number`. The final array has three dimensions, `forecastMonth`, `latitude` and `longitude`.

# In[35]:


seas5_anomalies_3m_202105_em = seas5_anomalies_3m_202105.mean('number')


# #### Convert precipitation rate to accumulation in mm
# 
# The last step before visualizing the 3-month seasonal forecast anomalies is again to convert the precipitation accumulations to total precipitation in mm. We repeat the same steps as previously:
# > 1. **Calculate number of days for each forecast month and add it as coordinate info** <br>
# > 2. **Name the 3-month rolling archives to indicate over which months the average was built** <br>
# > 3. **Convert the precipitation accumulations based on the number of days** <br>
# > 4. **Add updated attributes to the data array** <br>

# In[36]:


# Calculate number of days for each forecast month and add it as coordinate information to the data array
vt = [ pd.to_datetime(seas5_anomalies_3m_202105_em.time.values) + relativedelta(months=fcmonth-1) for fcmonth in seas5_anomalies_3m_202105_em.forecastMonth]
vts = [[thisvt+relativedelta(months=-mm) for mm in range(3)] for thisvt in vt]
numdays = [np.sum([monthrange(dd.year,dd.month)[1] for dd in d3]) for d3 in vts]
seas5_anomalies_3m_202105_em = seas5_anomalies_3m_202105_em.assign_coords(numdays=('forecastMonth',numdays))

# Define names for the 3-month rolling archives, that give an indication over which months the average was built
vts_names = ['{}{}{} {}'.format(d3[2].strftime('%b')[0],d3[1].strftime('%b')[0],d3[0].strftime('%b')[0], d3[0].strftime('%Y'))  for d3 in vts]
seas5_anomalies_3m_202105_em = seas5_anomalies_3m_202105_em.assign_coords(valid_time=('forecastMonth',vts_names))
seas5_anomalies_3m_202105_em

# Convert the precipitation accumulations based on the number of days
seas5_anomalies_3m_202105_em_tp = seas5_anomalies_3m_202105_em * seas5_anomalies_3m_202105_em.numdays * 24 * 60 * 60 * 1000

# Add updated attributes
seas5_anomalies_3m_202105_em_tp.attrs['units'] = 'mm'
seas5_anomalies_3m_202105_em_tp.attrs['long_name'] = 'SEAS5 3-monthly total precipitation ensemble mean anomaly for 6 lead-time months, start date in May 2021.'


# #### Visualise 3-monthly total precipitation ensemble mean anomalies
# 
# Let us now visualise the 3-month total precipitation ensemble mean anomaly in an interactive plot. Widgets (`ipywidgets`) allow us to add interactive features to Jupyter notebooks. We will use these to add a dropdown menu that offers the option to choose the 3-monthly periods over which the anomalies are averaged. Given that we have 6 lead-time months, we end up with only 4 possible 3-monthly periods (the last 2 months are insufficient to create a complete 3-month aggregation).

# In[37]:


vts_names


# In[38]:


dropdown_opts = [(vts_names[mm-1],mm) for mm in range(3,7)]
dropdown_opts


# In[39]:


tp_colors = [(153/255.,51/255.,0),(204/255.,136/255.,0),(1,213/255.,0),
             (1,238/255.,153/255.),(1,1,1),(204/255.,1,102/255.),
             (42/255.,1,0),(0,153/255.,51/255.),(0,102/255.,102/255.)]  
tp_levels = [-200,-100,-50,-10,10,50,100,200]



def plot_leadt(ds, leadt_end):
    array = ds.sel(forecastMonth=leadt_end)
    fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    im = plt.contourf(array.longitude, array.latitude, array, 
                      levels=tp_levels, colors=tp_colors, extend='both')
    ax.set_title(f'{array.long_name}' + os.linesep + f'{array.valid_time.values}', fontsize=16)
    ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--') 
    ax.coastlines(color='black')
    cbar = plt.colorbar(im,fraction=0.05, pad=0.04)
    cbar.set_label(array.units)
    
    return
    
dropdown_opts = [(vts_names[mm-1],mm) for mm in range(3,7)] 
a=widgets.interact(plot_leadt, ds=widgets.fixed(seas5_anomalies_3m_202105_em_tp), 
                   leadt_end=widgets.Dropdown(options=dropdown_opts,description='Valid Time:', 
                                              style={'description_width': 'initial'}))


# <div class="alert alert-block alert-success">
# <b>NOTE</b>: <br>
#     The plots shown here are similar but not identical to those on the <a href="https://climate.copernicus.eu/charts/c3s_seasonal/">Copernicus Climate Change Service (C3S) website</a>, specifically <a href="https://climate.copernicus.eu/charts/c3s_seasonal/c3s_seasonal_spatial_ecmf_rain_3m?facets=Parameters,precipitation%3BCentres,ECMWF&time=2021050100,744,2021060100&type=ensm&area=area08">here</a>. These plots are not identical however. The white areas in both plots, for example, have different meanings: those on the C3S website have a significance test to white-out non statistically-significant anomalies.</div>

# <hr>

# <p></p>
# <span style='float:right'><p style=\"text-align:right;\">This project is licensed under <a href="./LICENSE">APACHE License 2.0</a>. | <a href=\"https://github.com/ecmwf-projects/copernicus-training">View on GitHub</a></span>
