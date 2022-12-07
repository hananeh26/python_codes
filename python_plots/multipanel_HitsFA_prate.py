import os
#from netCDF4 import Dataset
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
#import netCDF4 as nc
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
#from netCDF4 import Dataset as netcdf_dataset
from matplotlib import colorbar, colors
from operator import itemgetter
import matplotlib.colors as mcolors
import sys
import pygrib
import cf2cdm
import cfgrib
from matplotlib.colors import from_levels_and_colors,Normalize
from scipy.signal import savgol_filter
import csv
import pandas as pd
from datetime import datetime, timedelta
import re
import datetime
import glob,os
import itertools
from matplotlib.font_manager import FontProperties
from numpy import mean
import statistics
from scipy import stats

########################################################################################

# This code is plotting precipitation rate at genesis time for FA and Hits.

########################################################################################

names=[]
genesis=[]
not_genesis=[]
genesis_time=[]
em=[]
not_genesis_FT=[]

dirs=u'/work/noaa/aoml-hafs1/galaka/noscrub/B220/'
sdir=u'/work/noaa/aoml-hafs1/hjafary/Dorian/dorian05l.2019082306/'
dires2=u'/work/noaa/aoml-hafs1/galaka/FOR_HANA/bdeck_invests/'
dires3=u'/work/noaa/aoml-hafs1/galaka/FOR_HANA/adeck_invests/'
picture_dir=u'/work/noaa/aoml-hafs1/hjafary/MED_Refl/'
grib_files='/work/noaa/aoml-hafs1/hjafary/hwrf_core/'
hwrf_atcf=u'/work/noaa/aoml-hafs1/galaka/noscrub/B220/'

regex = r'invest+[0-99]+l'
for x in glob.glob('/work/noaa/aoml-hafs1/galaka/noscrub/B220/*.atcfunix'):
    if re.search(regex, x):
        names.append(os.path.basename(x))

names.sort()
#print(names)

for sfile in names:
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    A_file=pd.read_csv(dirs+sfile,delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=headers,index_col=False)

    FP=[]
    for i in range(len(A_file.Vmax)):
        if A_file.Vmax[i]>=  34:   
            FP.append(A_file.Forecast_Period[i])

################################# Identifying if wind speed is above 34 for at least 12 consecutive hours ################

    count=0
    sub=[]
    cnt=0
#    FJ = [0,3,3,6,6,9,12]
    for j in range(0,len(FP)):
        if j==0:               # Always grab the first element to compare with the second value. 
            em.append(FP[j])
            count +=1

        elif FP[j]==FP[j-1]:
#        em.remove(FJ[j])
            cnt +=1

        elif FP[j]-FP[j-1]==3:  # If the second value is 3 more than the prevoius then append to the list. 
            em.append(FP[j])
                   # Also add the number to the list even if it gets repeated. 
#            sub.append(FJ[j]-FJ[j-1]) # Add the subtraction result between (j)-(j-1) to the sub list
            count +=1

        else:
            if count <= 4:     # If the total count is equal or less than 4, clear the list and start again. 
                em.clear()
                em.append(FP[j])

            else:
                count=1
                break


############################## Categorize files into genesis and not genesis also append genesis time ##################


    if len(em)>=5:    # If the total len of em is greater than 5 --> [0,3,6,9,12], it meets genesis requirement. 
        genesis.append(sfile)
        genesis_time.append(em[0])
#        print(sfile,", 'Generated into TC' ")
#        print(genesis)
    else:
        not_genesis.append(sfile)
    em.clear()

############################## Find hit, misses, false alarms, and correct negatives ##################################

slist=os.listdir(sdir)
swrf_out=[]
hwrf_core=[]
adeck_list=[]

################## Grab all the directories with hwrfprs-core #########################

for (root, dirs, files) in os.walk(grib_files):
    for directory in dirs:
        if directory.endswith('L'):
            swrf_out.append(grib_files+directory+'/')

################# Read all the grib files in  direcotry ############################
swrf_out.sort()
for lstfiles in swrf_out:
    for (root, dirs, files) in os.walk(lstfiles):
        for filename in files:
            if filename.endswith('.grb2'): #and int(filename.split('.')[0][-3:-1])<=40:
                hwrf_core.append(filename)

Misses1=[]
Misses=[]
Misses_time=[]
Correct_negatives1=[]
Hit=0
False_alarm_cases=[]
False_alarm_time=[]
Hit_genesis1=[]
Hit_genesis_time1=[]
Correct_negatives=[]
Correct_negatives_t=[]
Correct_negatives_time=[]
regex = r'invest+[0-40]'
genesis.sort()

for k in range(len(genesis)):
    if re.search(regex,genesis[k]):
        Hit_genesis1.append(genesis[k])
        Hit_genesis_time1.append(genesis_time[k])
        Hit += 1
    else:
        False_alarm_cases.append(genesis[k])
        False_alarm_time.append(genesis_time[k])


for l in range(len(not_genesis)):
    if re.search(regex, not_genesis[l]):
        Misses1.append(not_genesis[l])
    else:
        Correct_negatives1.append(not_genesis[l])


################################# False Alarm ################################

False_alarm_cases.sort()
false_sid=[]
false_date=[]
false_fhr=[]
for at in range(len(False_alarm_cases)):
    false_date.append(str(False_alarm_cases[at].split('.')[1][:10])) #2019082306
    false_fhr.append(str(False_alarm_time[at])) #48
    false_sid.append(str(False_alarm_cases[at].split('.')[0][-3:-1]))

false_file = []
for sid,date,fhr in zip(false_sid,false_date,false_fhr):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.core\.0p015\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
                    print(f'Found this GRIB2 file --> {filename}')
                    false_file.append(filename)
        if not false_file:
            print(f'Did not find a GRIB2 file --> {filename}')


###################################################################

# This section finds storm centers.

falselat_e=[]
falselon_e=[]
for c in range(len(False_alarm_cases)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+False_alarm_cases[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(0,len(a_file.Forecast_Period)-1):
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(False_alarm_time[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            falselat_e.append(int(ext_lat[:-1])/10.0)
            falselon_e.append(-(int(ext_lon[:-1])/10.0)+360)


####################################################################################

# Define a box around the storm center and find the average length of lat and lon.

falselat_len=[]
falselon_len=[]
for k in range(len(false_file)):
    FHR1 = int(false_file[k].split('.')[-2][1:4])
    FHR_date1=int(false_file[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+false_sid[k]+'L/'+false_date[k]+'/'+false_file[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'shortName': 'prate','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e[k]-3 and item  <= falselat_e[k]+3: # find all data points within the limits. 
            min_lat.append(index)  # append the indexes for those data points.
    falselat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e[k])-3 and itm <= (falselon_e[k])+3:
            min_lon.append(idx)
    falselon_len.append(len(min_lon))


###################################################################################

# This section finds lat, lon , and precip data points within a storm centered 3 x 3. 

falsecount_file=[]
falselat_c=[]
falselon_c=[]
falselons=[]
falselats=[]
falserain=[]

for k in range(len(false_file)):
    FHR1 = int(false_file[k].split('.')[-2][1:4])
    FHR_date1=int(false_file[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+false_sid[k]+'L/'+false_date[k]+'/'+false_file[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'shortName': 'prate','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e[k]-3 and item  <= falselat_e[k]+3:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e[k])-3 and itm <= (falselon_e[k])+3:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(falselat_len)

    # Compare the length of the lat for each invest file to the average length.

    if len(minlat) < a[0][0]-1:  # remove the case if it has 2 or more points less than the average.
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1) # append the next data point if it's less than the average
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)           # remove the last point if it has extra data points.    

    b=stats.mode(falselon_len)   

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1) 
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)

    # find lat, lon, precip. within the defined box.    

    falselons.append(lons1[minlon[:]])  
    falselats.append(lats1[minlat[:]])
    falserain.append((grbs.variables['prate'][minlat[:],minlon[:]])*3600)

    falsecount_file.append(false_file[k])
    falselat_c.append(falselat_e[k])
    falselon_c.append(falselon_e[k])


    

############################## Hits ################################

Hit_genesis=[]
Hit_genesis_time=[]

for w in range(len(Hit_genesis1)):
    if Hit_genesis1[w].startswith('invest06l.2019082500'):
        continue
    elif Hit_genesis1[w].startswith('invest07l.2018090212'):
        continue
    elif Hit_genesis1[w].startswith('invest10l.2018091212'):
        continue
    elif Hit_genesis1[w].startswith('invest15l.2018100812'):
        continue
    elif Hit_genesis1[w].startswith('invest15l.2018100818'):
        continue
    elif Hit_genesis1[w].startswith('invest15l.2018100900'):
        continue
    elif Hit_genesis1[w].startswith('invest15l.2018100906'):
        continue
    elif Hit_genesis1[w].startswith('invest17l.2020090612'):
        continue
    elif Hit_genesis1[w].startswith('invest18l.2020090700'):
        continue
    elif Hit_genesis1[w].startswith('invest27l.2020101800'):
        continue
    elif Hit_genesis1[w].startswith('invest25l.2020100118'):
        continue
    elif Hit_genesis1[w].startswith('invest12l.2019092012'):
        continue
    elif Hit_genesis1[w].startswith('invest12l.2019092112'):
        continue
    elif Hit_genesis1[w].startswith('invest27l.2020101806'):
        continue
    elif Hit_genesis1[w].startswith('invest30l.2020110818'):
        continue

    else:
        Hit_genesis.append(Hit_genesis1[w])
        Hit_genesis_time.append(Hit_genesis_time1[w])

all_sid=[]
all_date=[]
all_fhr=[]
FOUND_FILE=[]
for k in range(len(Hit_genesis)):    
    all_date.append(str(Hit_genesis[k].split('.')[1][:10])) #2019082306
    all_fhr.append(str(Hit_genesis_time[k])) #48
    all_sid.append(str(Hit_genesis[k].split('.')[0][-3:-1]))

FOUND_FILE = []
for sid,date,fhr in zip(all_sid,all_date,all_fhr):
#    print(f'Trying to match: sid={sid}, date={date}, fhr={fhr}')
#    ex= re.compile(r'[a-z]+'+str(storm_name)+str(storm_num)+'.'+str(FHR_date)+'.hwrfprs.core.0p015.f'+str(FHR)+'.grb2')
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.core\.0p015\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    FOUND_FILE.append(filename)
        if not FOUND_FILE:
            print('Did not find a GRIB2 file.')


###################################################################

# This section finds storm centers.

lat_e=[]
lon_e=[]
Hit_genesis.sort()
for c in range(len(Hit_genesis)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Hit_genesis[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(0,len(a_file.Forecast_Period)-1):
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(Hit_genesis_time[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            lat_e.append(int(ext_lat[:-1])/10.0)
            lon_e.append(-(int(ext_lon[:-1])/10.0)+360)

####################################################################################

# Define a box around the storm center and find the average length of lat and lon.

lats=[]
lons=[]
lat_len=[]
lon_len=[]
lat_c=[]
lon_c=[]
rain=[]
count_file=[]
FOUND_FILE.sort()


for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':{'shortName': 'prate','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)): 
        if item >= lat_e[k]-3 and item  <= lat_e[k]+3:  # find all data points within the limits. 
            min_lat.append(index)        # append the indexes for those data points.
    lat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-3 and itm <= (lon_e[k])+3:
            min_lon.append(idx)
    lon_len.append(len(min_lon))


###################################################################################

for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'shortName': 'prate','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= lat_e[k]-3 and item  <= lat_e[k]+3:
            min_lat.append(index)

    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-3 and itm <= (lon_e[k])+3:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(lat_len)

    if len(minlat) < a[0][0]-1:    # Compare the length of the lat for each invest file to the average length.
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(lon_len)

    if len(minlon) < b[0][0]-1:  # Compare the length of the lon for each invest file to the average length.
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)       

    lons.append(lons1[minlon[:]])
    lats.append(lats1[minlat[:]])
    rain.append((grbs.variables['prate'][minlat[:],minlon[:]])*3600)

    count_file.append(FOUND_FILE[k])
    lat_c.append(lat_e[k])
    lon_c.append(lon_e[k])

############################# Hit mean ################################

lon_e=np.mean(lon_c)-360     # storm centers average
lat_e=np.mean(lat_c)

rain_mean=np.nanmean(rain,axis=0)
y=np.nanmean(lats,axis=0)
x=np.nanmean(lons,axis=0)

############################ False Alarms mean ###########################

#falsecount_file.append(false_file[k])
false_lons=np.nanmean(falselons,axis=0)
false_lats=np.nanmean(falselats,axis=0)
false_rain=np.nanmean(falserain,axis=0)

false_latc=np.nanmean(falselat_c)
false_lonc=np.nanmean(falselon_c)-360

################################## Plotting ###########################

##################### Hits ###########################
#fig=plt.figure(figsize=(15,8))
fig, (ax,ax2)=plt.subplots((2,2),figsize=(20,10))
ax = fig.add_subplot(2,1,1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution="50m",linewidths=1)
ax.set_extent([lon_e-3,lon_e+3,lat_e-3,lat_e+3])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([lon_e-3,lon_e,lon_e+3])
gl.ylocator = mticker.FixedLocator([lat_e-3,lat_e,lat_e+3])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

clevs=[0,0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.25,1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
                7, 10., 15, 20,25,30,35,40,50,60,70]

cmap_data = np.array([[255, 255, 255],[127, 255, 0],[0,205 ,0],
                 [0,139 ,0], [20, 75 ,137],[27, 145, 255],
                 [0, 177, 235], [0, 240, 240], [139, 102, 210],
                 [144, 43, 235], [147, 0, 139], [135, 0, 0],
                 [203, 0, 0], [242, 62, 0], [255, 122, 0],
                 [204, 131, 0], [255, 215, 0], [255, 255, 0],
                 [255, 173, 184], [255, 0, 255], [139, 102, 210],
                 [144, 43, 235], [140, 0, 140], [17, 77, 140],
                 [255, 230, 220]],np.float32) / 255.0

cmap_rain, norm_rain = from_levels_and_colors(clevs, cmap_data, extend="max")
contour_c=ax.pcolormesh(x, y, rain_mean, cmap=cmap_rain,norm=norm_rain)

cb = fig.colorbar(contour_c, orientation="vertical",
                      pad=0.02, aspect=16, shrink=0.8)
cb.set_label('mm/hr',size=10)
cb.ax.tick_params(labelsize=10)
font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.set_title("Avg Precipitation Rate for Hits",fontsize=14,fontproperties=font)
ax.text(0.05,0.95, 'Total Cases:'+ str(len(count_file)),transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
ax.set_title('Core Data',fontsize='small',fontproperties=font ,loc='right')

##################  False Alarms ######################

ax2=fig.add_subplot(2,1,2)
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.coastlines(resolution="50m",linewidths=1)
ax2.set_extent([false_lonc-3,false_lonc+3,false_latc-3,false_latc+3])
#        print(lon_e[m])
#        print(lat_e[m])
gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([false_lonc-3,false_lonc,false_lonc+3])
gl.ylocator = mticker.FixedLocator([false_latc-3,false_latc,false_latc+3])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

clevs=[0,0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.25,1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
                7, 10., 15, 20,25,30,35,40,50,60,70]

cmap_data = np.array([[255, 255, 255],[127, 255, 0],[0,205 ,0],
                 [0,139 ,0], [20, 75 ,137],[27, 145, 255],
                 [0, 177, 235], [0, 240, 240], [139, 102, 210],
                 [144, 43, 235], [147, 0, 139], [135, 0, 0],
                 [203, 0, 0], [242, 62, 0], [255, 122, 0],
                 [204, 131, 0], [255, 215, 0], [255, 255, 0],
                 [255, 173, 184], [255, 0, 255], [139, 102, 210],
                 [144, 43, 235], [140, 0, 140], [17, 77, 140],
                 [255, 230, 220]],np.float32) / 255.0

cmap_rain, norm_rain = from_levels_and_colors(clevs, cmap_data, extend="max")
contour_c=ax2.pcolormesh(false_lons, false_lats, false_rain, cmap=cmap_rain,norm=norm_rain)
font = FontProperties()
font.set_family('serif')
font.set_name('Times')
ax2.text(0.05,0.95, 'Total Cases:'+ str(len(falsecount_file)),transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
ax2.set_title("Avg Precipitation Rate for False Alarms",fontsize='small',fontproperties=font )
ax2.set_title('Core Files',fontsize='small',fontproperties=font ,loc='right')

fig.savefig(picture_dir+'Avg_Precipitation_Rate_multipanel'+'.png')

