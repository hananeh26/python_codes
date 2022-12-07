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

# This code is plotting precipitation from 48 hours before up to hour 0 for FA and Hits.

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
picture_dir=u'/work/noaa/aoml-hafs1/hjafary/ws_picture/'
grib_files=u'/work/noaa/aoml-hafs1/hjafary/hwrf_global/'
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


Time=[-48,-24,0,24,48]
CN=[]
Misses1=[]
Misses=[]
Misses_time=[]
Correct_negatives1=[]
Hit=0
False_alarm_cases=[]
False_alarm_time=[]
False_alarm_cases1=[]
False_alarm_time1=[]
Hit_genesis1=[]
Hit_genesis_time1=[]
Hit_genesis=[]
Hit_genesis_time=[]
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
        False_alarm_cases1.append(genesis[k])
        False_alarm_time1.append(genesis_time[k])


for l in range(len(not_genesis)):
    if re.search(regex, not_genesis[l]):
        Misses1.append(not_genesis[l])
    else:
        Correct_negatives1.append(not_genesis[l])

########## seperate cases that exist 48 hours before #############

for l in range(len(Hit_genesis_time1)):
    if Hit_genesis_time1[l]-48 >= 0:
        Hit_genesis_time.append(Hit_genesis_time1[l])
        Hit_genesis.append(Hit_genesis1[l])

for k in range(len(False_alarm_time1)):
    if False_alarm_time1[k]-48 >= 0:
        False_alarm_time.append(False_alarm_time1[k])
        False_alarm_cases.append(False_alarm_cases1[k])


######################## False Alarm Genesis ########################

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
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.global\.0p25\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    false_file.append(filename)
        if not false_file:
            print(f'Did not find a GRIB2 file --> {filename}')

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

falselat_len=[]
falselon_len=[]
for k in range(len(false_file)):
    FHR1 = int(false_file[k].split('.')[-2][1:4])
    FHR_date1=int(false_file[k].split('.')[1][:10])

    grbs = xr.open_dataset(grib_files+false_sid[k]+'L/'+false_date[k]+'/'+false_file[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e[k]-5 and item  <= falselat_e[k]+5:
            min_lat.append(index)
    falselat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e[k])-5 and itm <= (falselon_e[k])+5:
            min_lon.append(idx)
    falselon_len.append(len(min_lon))


falsecount_file=[]
falselat_ct1=[]
falselon_ct1=[]
falselons1=[]
falselats1=[]

shearmag_fa=[]
u250_fa=[]
u850_fa=[]
v250_fa=[]
v850_fa=[]

new_fa=[]
new_fa_time=[]

for k in range(len(false_file)):
    FHR1 = int(false_file[k].split('.')[-2][1:4])
    FHR_date1=int(false_file[k].split('.')[1][:10])

    grbs = xr.open_dataset(grib_files+false_sid[k]+'L/'+false_date[k]+'/'+false_file[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    grbs2 = xr.open_dataset(grib_files+false_sid[k]+'L/'+false_date[k]+'/'+false_file[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'u','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e[k]-5 and item  <= falselat_e[k]+5:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e[k])-5 and itm <= (falselon_e[k])+5:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(falselat_len)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(falselon_len)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)


    falselons1.append(lons1[minlon[:]])
    falselats1.append(lats1[minlat[:]])
    falsecount_file.append(false_file[k])
    falselat_ct1.append(falselat_e[k])
    falselon_ct1.append(falselon_e[k])

###################################### Calculate vertical wind shear ###########################

    u850=(grbs2.variables['u'][6,minlat[:],minlon[:]])
    u250=(grbs2.variables['u'][20,minlat[:],minlon[:]])
    v850=(grbs.variables['v'][6,minlat[:],minlon[:]])
    v250=(grbs.variables['v'][20,minlat[:],minlon[:]])

    shearmag_fa.append(np.sqrt((u250-u850)**2+(v250-v850)**2))

    u250_fa.append(grbs2.variables['u'][20,minlat[:],minlon[:]])
    u850_fa.append(grbs2.variables['u'][6,minlat[:],minlon[:]])
    v250_fa.append(grbs.variables['v'][20,minlat[:],minlon[:]])
    v850_fa.append(grbs.variables['v'][6,minlat[:],minlon[:]])

############################################################################################

    new_fa.append(False_alarm_cases[k])
    new_fa_time.append(False_alarm_time[k])

falselon_ct=round((np.nanmean(falselon_ct1)-360),1)
falselat_ct=round((np.nanmean(falselat_ct1)),1)
falselat=np.nanmean(falselats1,axis=0)
falselon=np.nanmean(falselons1,axis=0)

fa_shearmag=np.nanmean(shearmag_fa,axis=0)
fa_u250=np.nanmean(u250_fa,axis=0)
fa_u850=np.nanmean(u850_fa,axis=0)
fa_v250=np.nanmean(v250_fa,axis=0)
fa_v850=np.nanmean(v850_fa,axis=0)

fa_ushear=fa_u250-fa_u850
fa_vshear=fa_v250-fa_v850
fa_speed=np.hypot(fa_ushear,fa_vshear)

fa_wsu=(fa_u250-fa_u850)/fa_speed
fa_wsv=(fa_v250-fa_v850)/fa_speed

            ############################ False Alarms Timing ###################

False_alarm_min24=[]
False_alarm_t_min24=[]
for k in range(len(new_fa_time)):
    if new_fa_time[k]-24 >= 0:  # 
        False_alarm_t_min24.append(new_fa_time[k]-24)
        False_alarm_min24.append(new_fa[k])

False_alarm_24=[]
False_alarm_t_24=[]
for k in range(len(new_fa_time)):
    if new_fa_time[k]+24 <= 126:  # 
        False_alarm_t_24.append(new_fa_time[k]+24)
        False_alarm_24.append(new_fa[k])

False_alarm_min48=[]
False_alarm_t_min48=[]
for k in range(len(new_fa_time)):
    if new_fa_time[k]-48 >= 0:  # 
        False_alarm_t_min48.append(new_fa_time[k]-48)
        False_alarm_min48.append(new_fa[k])

False_alarm_t_48=[]
False_alarm_48=[]
for k in range(len(new_fa_time)):
    if new_fa_time[k]+48 <= 126:  # 
        False_alarm_t_48.append(new_fa_time[k]+48)
        False_alarm_48.append(new_fa[k])

             ################### False Alarms -48hr ######################

FA=[]
False_alarm_min48.sort()
false_sid_min48=[]
false_date_min48=[]
false_fhr_min48=[]

falselon1_min48=[]
falselat1_min48=[]

false_min48=[]
falselat_ct1_min48=[]
falselon_ct1_min48=[]

for at in range(len(False_alarm_min48)):
    false_date_min48.append(str(False_alarm_min48[at].split('.')[1][:10])) #2019082306
    false_fhr_min48.append(str(False_alarm_t_min48[at])) #48
    false_sid_min48.append(str(False_alarm_min48[at].split('.')[0][-3:-1]))

false_file_min48 = []
for sid,date,fhr in zip(false_sid_min48,false_date_min48,false_fhr_min48):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.global\.0p25\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    false_file_min48.append(filename)
        if not false_file_min48:
            print(f'Did not find a GRIB2 file --> {filename}')

falselat_e_min48=[]
falselon_e_min48=[]
for c in range(len(False_alarm_min48)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+False_alarm_min48[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(0,len(a_file.Forecast_Period)-1):
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(False_alarm_t_min48[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            falselat_e_min48.append(int(ext_lat[:-1])/10.0)
            falselon_e_min48.append(-(int(ext_lon[:-1])/10.0)+360)

falselat_len_min48=[]
falselon_len_min48=[]
for k in range(len(false_file_min48)):
    FHR1 = int(false_file_min48[k].split('.')[-2][1:4])
    FHR_date1=int(false_file_min48[k].split('.')[1][:10])

    grbs = xr.open_dataset(grib_files+false_sid_min48[k]+'L/'+false_date_min48[k]+'/'+false_file_min48[k],engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e_min48[k]-5 and item  <= falselat_e_min48[k]+5:
            min_lat.append(index)
    falselat_len_min48.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e_min48[k])-5 and itm <= (falselon_e_min48[k])+5:
            min_lon.append(idx)
    falselon_len_min48.append(len(min_lon))

shearmag_famin48=[]
u250_famin48=[]
u850_famin48=[]
v250_famin48=[]
v850_famin48=[]

for k in range(len(false_file_min48)):
    FHR1 = int(false_file_min48[k].split('.')[-2][1:4])
    FHR_date1=int(false_file_min48[k].split('.')[1][:10])

    grbs = xr.open_dataset(grib_files+false_sid_min48[k]+'L/'+false_date_min48[k]+'/'+false_file_min48[k],engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    grbs2 = xr.open_dataset(grib_files+false_sid_min48[k]+'L/'+false_date_min48[k]+'/'+false_file_min48[k],engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'u','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e_min48[k]-5 and item  <= falselat_e_min48[k]+5:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e_min48[k])-5 and itm <= (falselon_e_min48[k])+5:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(falselat_len_min48)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(falselon_len_min48)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)

    falselon1_min48.append(lons1[minlon[:]])
    falselat1_min48.append(lats1[minlat[:]])
    false_min48.append(false_file_min48[k])
    falselat_ct1_min48.append(falselat_e_min48[k])
    falselon_ct1_min48.append(falselon_e_min48[k])

############################ calculate shear magnitude #####################

    u850=(grbs2.variables['u'][6,minlat[:],minlon[:]])
    u250=(grbs2.variables['u'][20,minlat[:],minlon[:]])
    v850=(grbs.variables['v'][6,minlat[:],minlon[:]])
    v250=(grbs.variables['v'][20,minlat[:],minlon[:]])

    shearmag_famin48.append(np.sqrt((u250-u850)**2+(v250-v850)**2))

########################## calculate shear direction ######################

    u250_famin48.append(grbs2.variables['u'][20,minlat[:],minlon[:]])
    u850_famin48.append(grbs2.variables['u'][6,minlat[:],minlon[:]])
    v250_famin48.append(grbs.variables['v'][20,minlat[:],minlon[:]])
    v850_famin48.append(grbs.variables['v'][6,minlat[:],minlon[:]])

famin48_u250=np.nanmean(u250_famin48,axis=0)
famin48_u850=np.nanmean(u850_famin48,axis=0)
famin48_v250=np.nanmean(v250_famin48,axis=0)
famin48_v850=np.nanmean(v850_famin48,axis=0)

famin48_ushear=famin48_u250-famin48_u850
famin48_vshear=famin48_v250-famin48_v850
famin48_speed=np.hypot(famin48_ushear,famin48_vshear)

famin48_wsu=(famin48_u250-famin48_u850)/famin48_speed
famin48_wsv=(famin48_v250-famin48_v850)/famin48_speed




falselon_ct_min48=round((np.nanmean(falselon_ct1_min48)-360),1)
falselat_ct_min48=round((np.nanmean(falselat_ct1_min48)),1)
falselat_min48=np.nanmean(falselat1_min48,axis=0)
falselon_min48=np.nanmean(falselon1_min48,axis=0)

famin48_shearmag=np.nanmean(shearmag_famin48,axis=0)

          ################### False Alarms -24hr ######################

False_alarm_min24.sort()
false_sid_min24=[]
false_date_min24=[]
false_fhr_min24=[]

falselon1_min24=[]
falselat1_min24=[]

false_min24=[]
falselat_ct1_min24=[]
falselon_ct1_min24=[]

for at in range(len(False_alarm_min24)):
    false_date_min24.append(str(False_alarm_min24[at].split('.')[1][:10])) #2019082306
    false_fhr_min24.append(str(False_alarm_t_min24[at])) #48
    false_sid_min24.append(str(False_alarm_min24[at].split('.')[0][-3:-1]))

false_file_min24 = []
for sid,date,fhr in zip(false_sid_min24,false_date_min24,false_fhr_min24):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.global\.0p25\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    false_file_min24.append(filename)
        if not false_file_min24:
            print(f'Did not find a GRIB2 file --> {filename}')

falselat_e_min24=[]
falselon_e_min24=[]
for c in range(len(False_alarm_min24)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+False_alarm_min24[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(0,len(a_file.Forecast_Period)-1):
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(False_alarm_t_min24[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            falselat_e_min24.append(int(ext_lat[:-1])/10.0)
            falselon_e_min24.append(-(int(ext_lon[:-1])/10.0)+360)

falselat_len_min24=[]
falselon_len_min24=[]
for k in range(len(false_file_min24)):
    FHR1 = int(false_file_min24[k].split('.')[-2][1:4])
    FHR_date1=int(false_file_min24[k].split('.')[1][:10])

    grbs = xr.open_dataset(grib_files+false_sid_min24[k]+'L/'+false_date_min24[k]+'/'+false_file_min24[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e_min24[k]-5 and item  <= falselat_e_min24[k]+5:
            min_lat.append(index)
    falselat_len_min24.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e_min24[k])-5 and itm <= (falselon_e_min24[k])+5:
            min_lon.append(idx)
    falselon_len_min24.append(len(min_lon))

shearmag_famin24=[]
u250_famin24=[]
u850_famin24=[]
v250_famin24=[]
v850_famin24=[]

for k in range(len(false_file_min24)):
    FHR1 = int(false_file_min24[k].split('.')[-2][1:4])
    FHR_date1=int(false_file_min24[k].split('.')[1][:10])

    grbs = xr.open_dataset(grib_files+false_sid_min24[k]+'L/'+false_date_min24[k]+'/'+false_file_min24[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    grbs2 = xr.open_dataset(grib_files+false_sid_min24[k]+'L/'+false_date_min24[k]+'/'+false_file_min24[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'u','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e_min24[k]-5 and item  <= falselat_e_min24[k]+5:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e_min24[k])-5 and itm <= (falselon_e_min24[k])+5:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(falselat_len_min24)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(falselon_len_min24)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)

    falselon1_min24.append(lons1[minlon[:]])
    falselat1_min24.append(lats1[minlat[:]])
    false_min24.append(false_file_min24[k])
    falselat_ct1_min24.append(falselat_e_min24[k])
    falselon_ct1_min24.append(falselon_e_min24[k])

    u850=(grbs2.variables['u'][6,minlat[:],minlon[:]])
    u250=(grbs2.variables['u'][20,minlat[:],minlon[:]])
    v850=(grbs.variables['v'][6,minlat[:],minlon[:]])
    v250=(grbs.variables['v'][20,minlat[:],minlon[:]])

    shearmag_famin24.append(np.sqrt((u250-u850)**2+(v250-v850)**2))

    u250_famin24.append(grbs2.variables['u'][20,minlat[:],minlon[:]])
    u850_famin24.append(grbs2.variables['u'][6,minlat[:],minlon[:]])
    v250_famin24.append(grbs.variables['v'][20,minlat[:],minlon[:]])
    v850_famin24.append(grbs.variables['v'][6,minlat[:],minlon[:]])

famin24_shearmag=np.nanmean(shearmag_famin24,axis=0)
famin24_u250=np.nanmean(u250_famin24,axis=0)
famin24_u850=np.nanmean(u850_famin24,axis=0)
famin24_v250=np.nanmean(v250_famin24,axis=0)
famin24_v850=np.nanmean(v850_famin24,axis=0)

falselon_ct_min24=round((np.nanmean(falselon_ct1_min24)-360),1)
falselat_ct_min24=round((np.nanmean(falselat_ct1_min24)),1)
falselat_min24=np.nanmean(falselat1_min24,axis=0)
falselon_min24=np.nanmean(falselon1_min24,axis=0)

famin24_ushear=famin24_u250-famin24_u850
famin24_vshear=famin24_v250-famin24_v850
famin24_speed=np.hypot(famin24_ushear,famin24_vshear)

famin24_wsu=(famin24_u250-famin24_u850)/famin24_speed
famin24_wsv=(famin24_v250-famin24_v850)/famin24_speed
 
############################## Hits ################################

Hit_genesis1=[]
Hit_genesis_time1=[]

for w in range(len(Hit_genesis)):
    if Hit_genesis[w].startswith('invest06l.2019082500'):
        continue
    elif Hit_genesis[w].startswith('invest07l.2018090212'):
        continue
    elif Hit_genesis[w].startswith('invest10l.2018091212'):
        continue
    elif Hit_genesis[w].startswith('invest15l.2018100812'):
        continue
    elif Hit_genesis[w].startswith('invest15l.2018100818'):
        continue
    elif Hit_genesis[w].startswith('invest15l.2018100900'):
        continue
    elif Hit_genesis[w].startswith('invest15l.2018100906'):
        continue
    elif Hit_genesis[w].startswith('invest17l.2020090612'):
        continue
    elif Hit_genesis[w].startswith('invest18l.2020090700'):
        continue
    elif Hit_genesis[w].startswith('invest27l.2020101800'):
        continue
    elif Hit_genesis[w].startswith('invest25l.2020100118'):
        continue
    elif Hit_genesis[w].startswith('invest12l.2019092012'):
        continue
    elif Hit_genesis[w].startswith('invest12l.2019092112'):
        continue
    elif Hit_genesis[w].startswith('invest27l.2020101806'):
        continue
    elif Hit_genesis[w].startswith('invest30l.2020110818'):
        continue

    else:
        Hit_genesis1.append(Hit_genesis[w])
        Hit_genesis_time1.append(Hit_genesis_time[w])



#####################################################################################

              ############################ Hit Genesis ###########################                         

all_sid=[]
all_date=[]
all_fhr=[]
FOUND_FILE=[]
for k in range(len(Hit_genesis1)):
    all_date.append(str(Hit_genesis1[k].split('.')[1][:10])) #2019082306
    all_fhr.append(str(Hit_genesis_time1[k])) #48
    all_sid.append(str(Hit_genesis1[k].split('.')[0][-3:-1]))

FOUND_FILE = []
for sid,date,fhr in zip(all_sid,all_date,all_fhr):
#    print(f'Trying to match: sid={sid}, date={date}, fhr={fhr}')
#    ex= re.compile(r'[a-z]+'+str(storm_name)+str(storm_num)+'.'+str(FHR_date)+'.hwrfprs.core.0p015.f'+str(FHR)+'.grb2')
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.global\.0p25\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    FOUND_FILE.append(filename)
        if not FOUND_FILE:
            print('Did not find a GRIB2 file.')


lat_e=[]
lon_e=[]
Hit_genesis1.sort()
for c in range(len(Hit_genesis1)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Hit_genesis1[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(0,len(a_file.Forecast_Period)-1):
#        if a_file.Forecast_Period[j+1] == a_file.Forecast_Period[j]:
#            continue
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(Hit_genesis_time1[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            lat_e.append(int(ext_lat[:-1])/10.0)
            lon_e.append(-(int(ext_lon[:-1])/10.0)+360)

lats=[]
lons=[]
lat_len=[]
lon_len=[]
lat_c=[]
lon_c=[]
rain=[]
count_file=[]
FOUND_FILE.sort()
#for plot_list in new_list:
for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], \
      engine='cfgrib',backend_kwargs={'filter_by_keys': \
        {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= lat_e[k]-5 and item  <= lat_e[k]+5:
            min_lat.append(index)
    lat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-5 and itm <= (lon_e[k])+5:
            min_lon.append(idx)
    lon_len.append(len(min_lon))

shearmag_hit=[]
u250_hit=[]
u850_hit=[]
v250_hit=[]
v850_hit=[]
new_hitgenesis=[]
new_hitgenesis_time=[]

for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    grbs2 = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                {'typeOfLevel': 'isobaricInhPa','shortName': 'u','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= lat_e[k]-5 and item  <= lat_e[k]+5:
            min_lat.append(index)

    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-5 and itm <= (lon_e[k])+5:
            min_lon.append(idx)
    minlat=min_lat
    minlon=min_lon
    a=stats.mode(lat_len)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(lon_len)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)

    lons.append(lons1[minlon[:]])
    lats.append(lats1[minlat[:]])
    count_file.append(FOUND_FILE[k])
    lat_c.append(lat_e[k])
    lon_c.append(lon_e[k])

    u850=(grbs2.variables['u'][6,minlat[:],minlon[:]])
    u250=(grbs2.variables['u'][20,minlat[:],minlon[:]])
    v850=(grbs.variables['v'][6,minlat[:],minlon[:]])
    v250=(grbs.variables['v'][20,minlat[:],minlon[:]])

    shearmag_hit.append(np.sqrt((u250-u850)**2+(v250-v850)**2))

    u250_hit.append(grbs2.variables['u'][20,minlat[:],minlon[:]])
    u850_hit.append(grbs2.variables['u'][6,minlat[:],minlon[:]])
    v250_hit.append(grbs.variables['v'][20,minlat[:],minlon[:]])
    v850_hit.append(grbs.variables['v'][6,minlat[:],minlon[:]])

    new_hitgenesis.append(Hit_genesis1[k])
    new_hitgenesis_time.append(Hit_genesis_time1[k])

hitlon_c=round((np.nanmean(lon_c)-360),1)
hitlat_c=round((np.nanmean(lat_c)),1)
hitlat=np.nanmean(lats,axis=0)
hitlon=np.nanmean(lons,axis=0)

hit_shearmag=np.nanmean(shearmag_hit,axis=0)

hit_u250=np.nanmean(u250_hit,axis=0)
hit_u850=np.nanmean(u850_hit,axis=0)
hit_v250=np.nanmean(v250_hit,axis=0)
hit_v850=np.nanmean(v850_hit,axis=0)

hit_ushear=hit_u250-hit_u850
hit_vshear=hit_v250-hit_v850
hit_speed=np.hypot(hit_ushear,hit_vshear)

hit_wsu=(hit_u250-hit_u850)/hit_speed
hit_wsv=(hit_v250-hit_v850)/hit_speed


########################### Hit Timing ##########################

Hit_genesis_24=[]
Hit_genesis_t_24=[]
for k in range(len(new_hitgenesis_time)):
    if new_hitgenesis_time[k]+24 <= 126:  # 
        Hit_genesis_t_24.append(new_hitgenesis_time[k]+24)
        Hit_genesis_24.append(new_hitgenesis[k])


Hit_genesis_min24=[]
Hit_genesis_t_min24=[]
for k in range(len(new_hitgenesis_time)):
    if new_hitgenesis_time[k]-24 >= 0:  # 
        Hit_genesis_t_min24.append(new_hitgenesis_time[k]-24)
        Hit_genesis_min24.append(new_hitgenesis[k])

Hit_genesis_48=[]
Hit_genesis_t_48=[]
for k in range(len(new_hitgenesis_time)):
    if new_hitgenesis_time[k]+48 <= 126:  # 
        Hit_genesis_t_48.append(new_hitgenesis_time[k]+48)
        Hit_genesis_48.append(new_hitgenesis[k])

Hit_genesis_min48=[]
Hit_genesis_t_min48=[]
for k in range(len(new_hitgenesis_time)):
    if new_hitgenesis_time[k]-48 >= 0:  # 
        Hit_genesis_t_min48.append(new_hitgenesis_time[k]-48)
        Hit_genesis_min48.append(new_hitgenesis[k])

               ################### Hit -48 #######################

H=[]
Hit_genesis_min48.sort()
Hit_sid_min48=[]
Hit_date_min48=[]
Hit_fhr_min48=[]

hitlon1_min48=[]
hitlat1_min48=[]

hit_min48=[]
hitlat_ct1_min48=[]
hitlon_ct1_min48=[]

for at in range(len(Hit_genesis_min48)):
    Hit_date_min48.append(str(Hit_genesis_min48[at].split('.')[1][:10])) #2019082306
    Hit_fhr_min48.append(str(Hit_genesis_t_min48[at])) #48
    Hit_sid_min48.append(str(Hit_genesis_min48[at].split('.')[0][-3:-1]))

Hit_file_min48 = []
for sid,date,fhr in zip(Hit_sid_min48,Hit_date_min48,Hit_fhr_min48):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.global\.0p25\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    Hit_file_min48.append(filename)
        if not Hit_file_min48:
            print(f'Did not find a GRIB2 file --> {filename}')

Hitlat_e_min48=[]
Hitlon_e_min48=[]
for c in range(len(Hit_genesis_min48)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Hit_genesis_min48[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(1,len(a_file.Forecast_Period)):
        if a_file.Forecast_Period[j] - a_file.Forecast_Period[j-1]==3 and \
               a_file.Forecast_Period[j-1] == int(str(Hit_genesis_t_min48[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            Hitlat_e_min48.append(int(ext_lat[:-1])/10.0)
            Hitlon_e_min48.append(-(int(ext_lon[:-1])/10.0)+360)

Hitlat_len_min48=[]
Hitlon_len_min48=[]

for k in range(len(Hit_file_min48)):
    FHR1 = int(Hit_file_min48[k].split('.')[-2][1:4])
    FHR_date1=int(Hit_file_min48[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+Hit_sid_min48[k]+'L/'+Hit_date_min48[k]+'/'+Hit_file_min48[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= Hitlat_e_min48[k]-5 and item  <= Hitlat_e_min48[k]+5:
            min_lat.append(index)
    Hitlat_len_min48.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (Hitlon_e_min48[k])-5 and itm <= (Hitlon_e_min48[k])+5:
            min_lon.append(idx)
    Hitlon_len_min48.append(len(min_lon))

shearmag_hitmin48=[]
u250_hitmin48=[]
u850_hitmin48=[]
v250_hitmin48=[]
v850_hitmin48=[]

for k in range(len(Hit_file_min48)):
    FHR1 = int(Hit_file_min48[k].split('.')[-2][1:4])
    FHR_date1=int(Hit_file_min48[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+Hit_sid_min48[k]+'L/'+Hit_date_min48[k]+'/'+Hit_file_min48[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    grbs2 = xr.open_dataset(grib_files+Hit_sid_min48[k]+'L/'+Hit_date_min48[k]+'/'+Hit_file_min48[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'u','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= Hitlat_e_min48[k]-5 and item  <= Hitlat_e_min48[k]+5:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (Hitlon_e_min48[k])-5 and itm <= (Hitlon_e_min48[k])+5:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(Hitlat_len_min48)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(Hitlon_len_min48)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)

    u850=(grbs2.variables['u'][6,minlat[:],minlon[:]])
    u250=(grbs2.variables['u'][20,minlat[:],minlon[:]])
    v850=(grbs.variables['v'][6,minlat[:],minlon[:]])
    v250=(grbs.variables['v'][20,minlat[:],minlon[:]])

    shearmag_hitmin48.append(np.sqrt((u250-u850)**2+(v250-v850)**2))

    u250_hitmin48.append(grbs2.variables['u'][20,minlat[:],minlon[:]])
    u850_hitmin48.append(grbs2.variables['u'][6,minlat[:],minlon[:]])
    v250_hitmin48.append(grbs.variables['v'][20,minlat[:],minlon[:]])
    v850_hitmin48.append(grbs.variables['v'][6,minlat[:],minlon[:]])

    hitlon1_min48.append(lons1[minlon[:]])
    hitlat1_min48.append(lats1[minlat[:]])
    hit_min48.append(Hit_file_min48[k])
    hitlat_ct1_min48.append(Hitlat_e_min48[k])
    hitlon_ct1_min48.append(Hitlon_e_min48[k])

hitlon_ct_min48=round((np.nanmean(hitlon_ct1_min48)-360),1)
hitlat_ct_min48=round((np.nanmean(hitlat_ct1_min48)),1)
hitlat_min48=np.nanmean(hitlat1_min48,axis=0)
hitlon_min48=np.nanmean(hitlon1_min48,axis=0)

hitmin48_shearmag=np.nanmean(shearmag_hitmin48,axis=0)
hitmin48_u250=np.nanmean(u250_hitmin48,axis=0)
hitmin48_u850=np.nanmean(u850_hitmin48,axis=0)
hitmin48_v250=np.nanmean(v250_hitmin48,axis=0)
hitmin48_v850=np.nanmean(v850_hitmin48,axis=0)

hitmin48_ushear=hitmin48_u250-hitmin48_u850
hitmin48_vshear=hitmin48_v250-hitmin48_v850
hitmin48_speed=np.hypot(hitmin48_ushear,hitmin48_vshear)

hitmin48_wsu=(hitmin48_u250-hitmin48_u850)/hitmin48_speed
hitmin48_wsv=(hitmin48_v250-hitmin48_v850)/hitmin48_speed

               ################### Hit -24 #######################

Hit_genesis_min24.sort()
Hit_sid_min24=[]
Hit_date_min24=[]
Hit_fhr_min24=[]

hitlon1_min24=[]
hitlat1_min24=[]

hit_min24=[]
hitlat_ct1_min24=[]
hitlon_ct1_min24=[]

for at in range(len(Hit_genesis_min24)):
    Hit_date_min24.append(str(Hit_genesis_min24[at].split('.')[1][:10])) #2019082306
    Hit_fhr_min24.append(str(Hit_genesis_t_min24[at])) #48
    Hit_sid_min24.append(str(Hit_genesis_min24[at].split('.')[0][-3:-1]))


Hit_file_min24 = []
for sid,date,fhr in zip(Hit_sid_min24,Hit_date_min24,Hit_fhr_min24):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.global\.0p25\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
#            adeck_list.append(kj)
#                    print(f'Found this GRIB2 file --> {filename}')
                    Hit_file_min24.append(filename)
        if not Hit_file_min24:
            print(f'Did not find a GRIB2 file --> {filename}')

Hitlat_e_min24=[]
Hitlon_e_min24=[]
for c in range(len(Hit_genesis_min24)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Hit_genesis_min24[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(0,len(a_file.Forecast_Period)-1):
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(Hit_genesis_t_min24[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            Hitlat_e_min24.append(int(ext_lat[:-1])/10.0)
            Hitlon_e_min24.append(-(int(ext_lon[:-1])/10.0)+360)

Hitlat_len_min24=[]
Hitlon_len_min24=[]
for k in range(len(Hit_file_min24)):
    FHR1 = int(Hit_file_min24[k].split('.')[-2][1:4])
    FHR_date1=int(Hit_file_min24[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+Hit_sid_min24[k]+'L/'+Hit_date_min24[k]+'/'+Hit_file_min24[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= Hitlat_e_min24[k]-5 and item  <= Hitlat_e_min24[k]+5:
            min_lat.append(index)
    Hitlat_len_min24.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (Hitlon_e_min24[k])-5 and itm <= (Hitlon_e_min24[k])+5:
            min_lon.append(idx)
    Hitlon_len_min24.append(len(min_lon))

shearmag_hitmin24=[]
u250_hitmin24=[]
u850_hitmin24=[]
v250_hitmin24=[]
v850_hitmin24=[]

for k in range(len(Hit_file_min24)):
    FHR1 = int(Hit_file_min24[k].split('.')[-2][1:4])
    FHR_date1=int(Hit_file_min24[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+Hit_sid_min24[k]+'L/'+Hit_date_min24[k]+'/'+Hit_file_min24[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                {'typeOfLevel': 'isobaricInhPa','shortName': 'v','stepRange':str(FHR1)},'indexpath':''})

    grbs2 = xr.open_dataset(grib_files+Hit_sid_min24[k]+'L/'+Hit_date_min24[k]+'/'+Hit_file_min24[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                {'typeOfLevel': 'isobaricInhPa','shortName': 'u','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= Hitlat_e_min24[k]-5 and item  <= Hitlat_e_min24[k]+5:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (Hitlon_e_min24[k])-5 and itm <= (Hitlon_e_min24[k])+5:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(Hitlat_len_min24)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(Hitlon_len_min24)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)
  
    u850=(grbs2.variables['u'][6,minlat[:],minlon[:]])
    u250=(grbs2.variables['u'][20,minlat[:],minlon[:]])
    v850=(grbs.variables['v'][6,minlat[:],minlon[:]])
    v250=(grbs.variables['v'][20,minlat[:],minlon[:]])

    shearmag_hitmin24.append(np.sqrt((u250-u850)**2+(v250-v850)**2))

    u250_hitmin24.append(grbs2.variables['u'][20,minlat[:],minlon[:]])
    u850_hitmin24.append(grbs2.variables['u'][6,minlat[:],minlon[:]])
    v250_hitmin24.append(grbs.variables['v'][20,minlat[:],minlon[:]])
    v850_hitmin24.append(grbs.variables['v'][6,minlat[:],minlon[:]])

    hitlon1_min24.append(lons1[minlon[:]])
    hitlat1_min24.append(lats1[minlat[:]])
    hit_min24.append(Hit_file_min24[k])
    hitlat_ct1_min24.append(Hitlat_e_min24[k])
    hitlon_ct1_min24.append(Hitlon_e_min24[k])

hitlon_ct_min24=round((np.nanmean(hitlon_ct1_min24)-360),1)
hitlat_ct_min24=round((np.nanmean(hitlat_ct1_min24)),1)
hitlat_min24=np.nanmean(hitlat1_min24,axis=0)
hitlon_min24=np.nanmean(hitlon1_min24,axis=0)

hitmin24_shearmag=np.nanmean(shearmag_hitmin24,axis=0)
hitmin24_u250=np.nanmean(u250_hitmin24,axis=0)
hitmin24_u850=np.nanmean(u850_hitmin24,axis=0)
hitmin24_v250=np.nanmean(v250_hitmin24,axis=0)
hitmin24_v850=np.nanmean(v850_hitmin24,axis=0)

hitmin24_ushear=hitmin24_u250-hitmin24_u850
hitmin24_vshear=hitmin24_v250-hitmin24_v850
hitmin24_speed=np.hypot(hitmin24_ushear,hitmin24_vshear)

hitmin24_wsu=(hitmin24_u250-hitmin24_u850)/hitmin24_speed
hitmin24_wsv=(hitmin24_v250-hitmin24_v850)/hitmin24_speed

################################## Plotting ###########################

     ####################### -48hr #######################


#fig.suptitle("Mean of 500-850mb Vertical Wind Shear 2017-2020 (-48hr to Genesis)",fontsize=20)
fig1 = plt.figure(figsize=(10,7))
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.coastlines(resolution="50m",linewidths=1)
ax1.set_extent([hitlon_ct_min48-5.01,hitlon_ct_min48+5,hitlat_ct_min48-5,hitlat_ct_min48+5.01])
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([hitlon_ct_min48-5,hitlon_ct_min48,hitlon_ct_min48+5])
gl.ylocator = mticker.FixedLocator([hitlat_ct_min48-5,hitlat_ct_min48,hitlat_ct_min48+5])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

clevs=np.arange(0, 16, 1.)

cmap_data = np.array([[255, 255, 255],[0,205 ,0],
                 [0,139 ,0], [20, 75 ,137],[27, 145, 255],
                 [0, 177, 235], [0, 240, 240], [139, 102, 210],
                 [144, 43, 235], [147, 0, 139], [135, 0, 0],
                 [242, 62, 0], [255, 122, 0],
                 [204, 131, 0], [255, 215, 0], [255, 255, 0]],np.float32) / 255.0


cmap_rain, norm_rain = from_levels_and_colors(clevs, cmap_data,
                                           extend="max")

#work from here
contour_c=ax1.pcolormesh(hitlon_min48, hitlat_min48, hitmin48_shearmag,cmap=cmap_rain,norm=norm_rain)

ax1.quiver(hitlon_min48[::3], hitlat_min48[::3], hitmin48_wsu[::3,::3], hitmin48_wsv[::3,::3], \
            headlength=4, headwidth=3, color= 'k',scale=None,label='300-700mb shear')

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.set_title("Hits -48hr",fontsize=12,fontproperties=font)
ax1.text(0.05,0.95, 'Total Cases:'+ str(len(hit_min48)),transform=ax1.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
cb = fig1.colorbar(contour_c, ax=[ax1], \
orientation="vertical",\
                      pad=0.02, aspect=10, shrink=0.8)
cb.set_label('m/s',size=10)
cb.ax.tick_params(labelsize=10)
fig1.savefig(picture_dir+'MidWS_min48hits'+'.png')



fig2 = plt.figure(figsize=(10,7))
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.set_extent([falselon_ct_min48-5.01,falselon_ct_min48+5,falselat_ct_min48-5,falselat_ct_min48+5.01])
gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([falselon_ct_min48-5,falselon_ct_min48,falselon_ct_min48+5])
gl.ylocator = mticker.FixedLocator([falselat_ct_min48-5,falselat_ct_min48,falselat_ct_min48+5])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

contour_c=ax2.pcolormesh(falselon_min48, falselat_min48, famin48_shearmag, cmap=cmap_rain,norm=norm_rain)
ax2.quiver(falselon_min48[::3], falselat_min48[::3], famin48_wsu[::3,::3], famin48_wsv[::3,::3], \
            headlength=4, headwidth=3, color= 'k',scale=None,label='300-700mb shear')

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax2.set_title("FA -48hr",fontsize=12,fontproperties=font)
ax2.text(0.05,0.95, 'Total Cases:'+ str(len(false_min48)),transform=ax2.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
cb = fig2.colorbar(contour_c, ax=[ax2], \
orientation="vertical",\
                      pad=0.02, aspect=10, shrink=0.8)
cb.set_label('m/s',size=10)
cb.ax.tick_params(labelsize=10)
fig2.savefig(picture_dir+'MidWS_min24FA'+'.png')



fig3 = plt.figure(figsize=(10,7))
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.set_extent([hitlon_ct_min24-5.01,hitlon_ct_min24+5,hitlat_ct_min24-5,hitlat_ct_min24+5.01])
gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([hitlon_ct_min24-5,hitlon_ct_min24,hitlon_ct_min24+5])
gl.ylocator = mticker.FixedLocator([hitlat_ct_min24-5,hitlat_ct_min24,hitlat_ct_min24+5])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}
contour_c=ax3.pcolormesh(hitlon_min24, hitlat_min24, hitmin24_shearmag,cmap=cmap_rain,norm=norm_rain)
ax3.quiver(hitlon_min24[::3], hitlat_min24[::3], hitmin24_wsu[::3,::3], hitmin24_wsv[::3,::3], \
            headlength=4, headwidth=3, color= 'k',scale=None,label='300-700mb shear')

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax3.set_title("Hits -24hr",fontsize=12,fontproperties=font)
ax3.text(0.05,0.95, 'Total Cases:'+ str(len(hit_min24)),transform=ax3.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
cb = fig3.colorbar(contour_c, ax=[ax3], \
orientation="vertical",\
                      pad=0.02, aspect=10, shrink=0.8)
cb.set_label('m/s',size=10)
cb.ax.tick_params(labelsize=10)
fig3.savefig(picture_dir+'MidWS_min24hits'+'.png')


fig4 = plt.figure(figsize=(10,7))
ax4 = plt.axes(projection=ccrs.PlateCarree())
ax4.coastlines(resolution="50m",linewidths=1)
ax4.set_extent([falselon_ct_min24-5.01,falselon_ct_min24+5,falselat_ct_min24-5,falselat_ct_min24+5.01])
gl = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([falselon_ct_min24-5,falselon_ct_min24,falselon_ct_min24+5])
gl.ylocator = mticker.FixedLocator([falselat_ct_min24-5,falselat_ct_min24,falselat_ct_min24+5])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

contour_c=ax4.pcolormesh(falselon_min24, falselat_min24, famin24_shearmag, cmap=cmap_rain,norm=norm_rain)
ax4.quiver(falselon_min24[::3], falselat_min24[::3], famin24_wsu[::3,::3], famin24_wsv[::3,::3], \
            headlength=4, headwidth=3, color= 'k',scale=None,label='300-700mb shear')

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax4.set_title("FA -24hr",fontsize=12,fontproperties=font)
ax4.text(0.05,0.95, 'Total Cases:'+ str(len(false_min24)),transform=ax4.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
cb = fig4.colorbar(contour_c, ax=[ax4], \
orientation="vertical",\
                      pad=0.02, aspect=10, shrink=0.8)
cb.set_label('m/s',size=10)
cb.ax.tick_params(labelsize=10)
fig4.savefig(picture_dir+'MidWS_min24FA'+'.png')


fig5 = plt.figure(figsize=(10,7))
ax5 = plt.axes(projection=ccrs.PlateCarree())
ax5.coastlines(resolution="50m",linewidths=1)
ax5.set_extent([hitlon_c-5.01,hitlon_c+5,hitlat_c-5,hitlat_c+5.01])
gl = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([hitlon_c-5,hitlon_c,hitlon_c+5])
gl.ylocator = mticker.FixedLocator([hitlat_c-5,hitlat_c,hitlat_c+5])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}
contour_c=ax5.pcolormesh(hitlon, hitlat, hit_shearmag, cmap=cmap_rain,norm=norm_rain)
ax5.quiver(hitlon[::3], hitlat[::3], hit_wsu[::3,::3], hit_wsv[::3,::3], \
            headlength=4, headwidth=3, color= 'k',scale=None,label='300-700mb shear')

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax5.set_title("Hits at Genesis",fontsize=12,fontproperties=font)
ax5.text(0.05,0.95, 'Total Cases:'+ str(len(count_file)),transform=ax5.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
cb = fig5.colorbar(contour_c, ax=[ax5], \
orientation="vertical",\
                      pad=0.02, aspect=10, shrink=0.8)
cb.set_label('m/s',size=10)
cb.ax.tick_params(labelsize=10)
fig5.savefig(picture_dir+'MidWS_minhits'+'.png')


fig6 = plt.figure(figsize=(10,7))
ax6 = plt.axes(projection=ccrs.PlateCarree())
ax6.coastlines(resolution="50m",linewidths=1)
ax6.set_extent([falselon_ct-5.01,falselon_ct+5,falselat_ct-5,falselat_ct+5.01])
#        print(lon_e[m])
#        print(lat_e[m])
gl = ax6.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([falselon_ct-5,falselon_ct,falselon_ct+5])
gl.ylocator = mticker.FixedLocator([falselat_ct-5,falselat_ct,falselat_ct+5])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

contour_c=ax6.pcolormesh(falselon, falselat,fa_shearmag, cmap=cmap_rain,norm=norm_rain)
ax6.quiver(falselon[::3], falselat[::3], fa_wsu[::3,::3], fa_wsv[::3,::3], \
            headlength=4, headwidth=3, color= 'k',scale=None,label='300-700mb shear')

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
ax6.text(0.05,0.95, 'Total Cases:'+ str(len(falsecount_file)),transform=ax6.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax6.set_title("FA at Genesis",fontsize=12,fontproperties=font )
cb = fig6.colorbar(contour_c, ax=[ax6], \
orientation="vertical",\
                      pad=0.02, aspect=10, shrink=0.8)
cb.set_label('m/s',size=10)
cb.ax.tick_params(labelsize=10)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
fig6.savefig(picture_dir+'DeepWS_minFA'+'.png')
                                               


