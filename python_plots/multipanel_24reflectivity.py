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

# This code is plotting composite reflectivity from 48 hours before up to hour 0 for all the cases.
# It is taking the median instead of average. 

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

temp_hit_time=[]
temp_hit=[]
temp_fa_cases=[]
temp_fa_time=[]

for k in range(len(genesis)):
    if re.search(regex,genesis[k]):
        temp_hit.append(genesis[k])
        temp_hit_time.append(genesis_time[k])
#        Hit_genesis1.append(genesis[k])
#        Hit_genesis_time1.append(genesis_time[k])
    else:
        temp_fa_cases.append(genesis[k])
        temp_fa_time.append(genesis_time[k])
#        False_alarm_cases.append(genesis[k])
#        False_alarm_time.append(genesis_time[k])


for l in range(len(not_genesis)):
    if re.search(regex, not_genesis[l]):
        Misses1.append(not_genesis[l])
    else:
        Correct_negatives1.append(not_genesis[l])

for k in range(len(temp_hit_time)):
    if temp_hit_time[k]+24 <= 126:  # 
        Hit_genesis_time1.append(temp_hit_time[k]+24)
        Hit_genesis1.append(temp_hit[k])

for k in range(len(temp_fa_time)):
    if temp_fa_time[k]+24 <= 126:  # 
        False_alarm_time.append(temp_fa_time[k]+24)
        False_alarm_cases.append(temp_fa_cases[k])

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

falselat_e=[]
falselon_e=[]
for c in range(len(False_alarm_cases)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+False_alarm_cases[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(1,len(a_file.Forecast_Period)):
        if a_file.Forecast_Period[j] - a_file.Forecast_Period[j-1]==3 and \
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
               {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= falselat_e[k]-3 and item  <= falselat_e[k]+3:
            min_lat.append(index)
    falselat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (falselon_e[k])-3 and itm <= (falselon_e[k])+3:
            min_lon.append(idx)
    falselon_len.append(len(min_lon))

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
               {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})
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


    falselons.append(lons1[minlon[:]])
    falselats.append(lats1[minlat[:]])
    falserain.append((grbs.variables['refd'][12,minlat[:],minlon[:]])) # 700mb

    falsecount_file.append(false_file[k])
    falselat_c.append(falselat_e[k])
    falselon_c.append(falselon_e[k])


    
################################# Misses ######################################

dirs=u'/work/noaa/aoml-hafs1/galaka/noscrub/B220/'
Misses_t=[]
temp_miss_t=[]
temp_miss=[]
for gfile in Misses1:
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    Gfile=pd.read_csv(dirs+gfile,delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=headers,index_col=False)
    for t in range(len(Gfile.Vmax)):
        if max(Gfile.Forecast_Period) == 126:
            if Gfile.Vmax[t]>= 34 or Gfile.Vmax[t]==max(Gfile.Vmax):
                temp_miss_t.append(Gfile.Forecast_Period[t])   
                temp_miss.append(gfile)     
#            Misses_t.append(Gfile.Forecast_Period[t])
#            Misses.append(gfile)
                break


for k in range(len(temp_miss_t)):
    if temp_miss_t[k]+24 <= 126:
        Misses_t.append(temp_miss_t[k]+24)
        Misses.append(temp_miss[k])


miss_sid=[]
miss_date=[]
miss_fhr=[]
for at in range(len(Misses)):
    miss_date.append(str(Misses[at].split('.')[1][:10])) #2019082306
    miss_fhr.append(str(Misses_t[at])) #48
    miss_sid.append(str(Misses[at].split('.')[0][-3:-1]))

miss_file1 = []
for sid,date,fhr in zip(miss_sid,miss_date,miss_fhr):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.core\.0p015\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
                    miss_file1.append(filename)
        if not miss_file1:
            print(f'Did not find a GRIB2 file --> {filename}')


misslat_e=[]
misslon_e=[]
for c in range(len(Misses)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Misses[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(1,len(a_file.Forecast_Period)):
#        if a_file.Forecast_Period[j+1] == a_file.Forecast_Period[j]:
#            continue
        if a_file.Forecast_Period[j] - a_file.Forecast_Period[j-1]==3 and \
               a_file.Forecast_Period[j] == int(str(Misses_t[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            misslat_e.append(int(ext_lat[:-1])/10.0)
            misslon_e.append(-(int(ext_lon[:-1])/10.0)+360)


misslat_len=[]
misslon_len=[]
miss_file1.sort()
for k in range(len(miss_file1)):
    FHR1 = int(miss_file1[k].split('.')[-2][1:4])
    FHR_date1=str(miss_file1[k].split('.')[1][:10])
    sid=str(miss_file1[k].split('.')[0][-3:-1])
    grbs = xr.open_dataset(grib_files+sid+'L/'+FHR_date1+'/'+miss_file1[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                   {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= misslat_e[k]-3 and item  <= misslat_e[k]+3:
            min_lat.append(index)
    if len(min_lat) != 0:
        misslat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (misslon_e[k])-3 and itm <= (misslon_e[k])+3:
            min_lon.append(idx)
    if len(min_lon) != 0:
        misslon_len.append(len(min_lon))

misscount_file=[]
misslat_c=[]
misslon_c=[]
misslons=[]
misslats=[]
missrain=[]

for k in range(len(miss_file1)):
    FHR1 = int(miss_file1[k].split('.')[-2][1:4])
    FHR_date1=str(miss_file1[k].split('.')[1][:10])
    sid=str(miss_file1[k].split('.')[0][-3:-1])
    grbs = xr.open_dataset(grib_files+sid+'L/'+FHR_date1+'/'+miss_file1[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                   {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= misslat_e[k]-3 and item  <= misslat_e[k]+3:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (misslon_e[k])-3 and itm <= (misslon_e[k])+3:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(misslat_len)
    b=stats.mode(misslon_len)

    temp_lons=(lons1[minlon[:]])
    temp_lats=(lats1[minlat[:]])

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0] and minlat[-1] < 600:
        minlat.append(minlat[-1]+1)
    elif len(minlat) < a[0][0] and minlat[-1] == 600:
        minlat.insert(0,minlat[0]-1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0] and minlon[-1] < 600:
        minlon.append(minlon[-1]+1)
    elif len(minlon) < b[0][0] and minlon[-1] == 600:
        minlon.insert(0,minlon[0]-1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)

    misslons.append(lons1[minlon[:]])
    misslats.append(lats1[minlat[:]])
    missrain.append((grbs.variables['refd'][12,minlat[:],minlon[:]]))

    misscount_file.append(miss_file1[k])
    misslat_c.append(misslat_e[k])
    misslon_c.append(misslon_e[k])


############################### Correct Negatives ############################
temp_cn1=[]
temp_cn_t1=[]
temp_cn=[]
temp_cn_t=[]
dirs=u'/work/noaa/aoml-hafs1/galaka/noscrub/B220/'
for gfile in Correct_negatives1:
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    Gfile=pd.read_csv(dirs+gfile,delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=headers,index_col=False)
# read the files in not_genesis list

    for t in range(len(Gfile.Vmax)):
        if max(Gfile.Forecast_Period) == 126:
            if Gfile.Vmax[t]>= 34 or Gfile.Vmax[t]==max(Gfile.Vmax):
                temp_cn_t1.append(Gfile.Forecast_Period[t])
                temp_cn1.append(gfile)
#            Correct_negatives_t.append(Gfile.Forecast_Period[t])
#            Correct_negatives.append(gfile)
                break

Correct_negatives_t=[]
Correct_negatives=[]
for k in range(len(temp_cn_t1)):
    if temp_cn_t1[k]+24 <= 126:
        Correct_negatives_t.append(temp_cn_t1[k]+24)
        Correct_negatives.append(temp_cn1[k])

cn_sid=[]
cn_date=[]
cn_fhr=[]
cn_file1=[]
cn_file=[]

for at in range(len(Correct_negatives)):
    cn_date.append(str(Correct_negatives[at].split('.')[1][:10])) #2019082306
    cn_fhr.append(str(Correct_negatives_t[at])) #48
    cn_sid.append(str(Correct_negatives[at].split('.')[0][-3:-1]))

for sid,date,fhr in zip(cn_sid,cn_date,cn_fhr):
    ex = re.compile(r'[a-z]+'+sid+'l'+'\.'+date+'\.hwrfprs\.core\.0p015\.f'+fhr.zfill(3)+'\.grb2')
    for lstfiles in swrf_out:
        for (root, dirs, files) in os.walk(lstfiles):
            for filename in files:
                if ex.match(filename):
                    cn_file1.append(filename)
            if not cn_file1:
                print(f'Did not find a GRIB2 file --> {filename}')


cnlat_e=[]
cnlon_e=[]
for c in range(len(Correct_negatives)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Correct_negatives[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(1,len(a_file.Forecast_Period)):
#        if a_file.Forecast_Period[j+1] == a_file.Forecast_Period[j]:
#            continue
        if a_file.Forecast_Period[j] - a_file.Forecast_Period[j-1]==3 and \
               a_file.Forecast_Period[j] == int(str(Correct_negatives_t[c]).zfill(3)):
            ext_lat=a_file.LatN[j]
            ext_lon=a_file.LonW[j]
            cnlat_e.append(int(ext_lat[:-1])/10.0)
            cnlon_e.append(-(int(ext_lon[:-1])/10.0)+360)

cnlat_len=[]
cnlon_len=[]

for k in range(len(cn_file1)):
    FHR1 = int(cn_file1[k].split('.')[-2][1:4])
    FHR_date1=str(cn_file1[k].split('.')[1][:10])
    sid=str(cn_file1[k].split('.')[0][-3:-1])
    grbs = xr.open_dataset(grib_files+sid+'L/'+FHR_date1+'/'+cn_file1[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                   {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])
#    print(lats1)
    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= cnlat_e[k]-3 and item  <= cnlat_e[k]+3:
            min_lat.append(index)
    cnlat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (cnlon_e[k])-3 and itm <= (cnlon_e[k])+3:
            min_lon.append(idx)
    cnlon_len.append(len(min_lon))

cnlons=[]
cnlats=[]
cnrain=[]
cncount_file=[]
cnlat_c=[]
cnlon_c=[]

for k in range(len(cn_file1)):
    FHR1 = int(cn_file1[k].split('.')[-2][1:4])
    FHR_date1=str(cn_file1[k].split('.')[1][:10])
    sid=str(cn_file1[k].split('.')[0][-3:-1])
    grbs = xr.open_dataset(grib_files+sid+'L/'+FHR_date1+'/'+cn_file1[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
                   {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= cnlat_e[k]-3 and item  <= cnlat_e[k]+3:
            min_lat.append(index)


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (cnlon_e[k])-3 and itm <= (cnlon_e[k])+3:
            min_lon.append(idx)

    minlat=min_lat
    minlon=min_lon
    a=stats.mode(cnlat_len)

    if len(minlat) < a[0][0]-1:
        continue
    elif len(minlat) < a[0][0]:
        minlat.append(minlat[-1]+1)
    elif len(minlat) > a[0][0]:
        minlat.pop(-1)

    b=stats.mode(cnlon_len)

    if len(minlon) < b[0][0]-1:
        continue
    elif len(minlon) < b[0][0]:
        minlon.append(minlon[-1]+1)
    elif len(minlon) > b[0][0]:
        minlon.pop(-1)
    
    cnlons.append(lons1[minlon[:]])
    cnlats.append(lats1[minlat[:]])
    cnrain.append((grbs.variables['refd'][12,minlat[:],minlon[:]]))

    cncount_file.append(cn_file1[k])
    cnlat_c.append(cnlat_e[k])
    cnlon_c.append(cnlon_e[k])

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


lat_e=[]
lon_e=[]
Hit_genesis.sort()
for c in range(len(Hit_genesis)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+Hit_genesis[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)
    for j in range(1,len(a_file.Forecast_Period)):
#        if a_file.Forecast_Period[j+1] == a_file.Forecast_Period[j]:
#            continue
        if a_file.Forecast_Period[j] - a_file.Forecast_Period[j-1]==3 and \
               a_file.Forecast_Period[j] == int(str(Hit_genesis_time[c]).zfill(3)):
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
    engine='cfgrib',backend_kwargs={'filter_by_keys':{'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])

    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= lat_e[k]-3 and item  <= lat_e[k]+3:
            min_lat.append(index)
    lat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-3 and itm <= (lon_e[k])+3:
            min_lon.append(idx)
    lon_len.append(len(min_lon))

for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'shortName': 'refd','stepRange':str(FHR1)},'indexpath':''})

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
    rain.append((grbs.variables['refd'][12,minlat[:],minlon[:]]))

    count_file.append(FOUND_FILE[k])
    lat_c.append(lat_e[k])
    lon_c.append(lon_e[k])
############################# Hit mean ################################

lon_e=np.median(lon_c)-360
lat_e=np.median(lat_c)

rain_mean=np.nanmedian(rain,axis=0)
y=np.nanmedian(lats,axis=0)
x=np.nanmedian(lons,axis=0)
############################ Correct Negatives mean #########################

#cncount_file.append(cn_file[k])
cn_lons=np.nanmedian(cnlons,axis=0)
cn_lats=np.nanmedian(cnlats,axis=0)
cn_rain=np.nanmedian(cnrain,axis=0)

cn_latc=np.nanmedian(cnlat_c)
cn_lonc=np.nanmedian(cnlon_c)-360

############################ Misses mean #################################

#misscount_file.append(miss_file[k])
miss_lons=np.nanmedian(misslons,axis=0)
miss_lats=np.nanmedian(misslats,axis=0)
miss_rain=np.nanmedian(missrain,axis=0)

miss_latc=np.nanmedian(misslat_c)
miss_lonc=np.nanmedian(misslon_c)-360

############################ False Alarms mean ###########################

#falsecount_file.append(false_file[k])
false_lons=np.nanmedian(falselons,axis=0)
false_lats=np.nanmedian(falselats,axis=0)
false_rain=np.nanmedian(falserain,axis=0)

false_latc=np.nanmedian(falselat_c)
false_lonc=np.nanmedian(falselon_c)-360

################################## Plotting ###########################

##################### Hits ###########################
fig=plt.figure(figsize=(14, 10))
fig.suptitle("700mb Median Refl +24hr",fontsize=18)
ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
ax1.coastlines(resolution="50m",linewidths=1)
ax1.set_extent([lon_e-3,lon_e+3,lat_e-3,lat_e+3])
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
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

dbz_levels = [1,5,9,13,17,21,25,29,33,37,41,45,49,53]

dbz_rgb = np.array([[4,233,231],
                    [1,159,244], [3,0,244],
                    [2,253,2], [1,197,1],
                    [0,142,0], [253,248,2],
                    [229,188,0], [253,149,0],
                    [253,0,0], [212,0,0],
                    [188,0,0],[248,0,253],
                    [152,84,198]], np.float32) / 255.0
dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb,\
                                           extend="max")
contour_c=ax1.pcolormesh(x, y, rain_mean,cmap=dbz_map,norm=dbz_norm)

font = FontProperties()
font.set_family('serif')
font.set_name('Times')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.set_title("700mb Median Refl for Hits",fontsize=14,fontproperties=font)
ax1.text(0.05,0.95, 'Total Cases:'+ str(len(count_file)),transform=ax1.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
#ax.title('Core Data',fontsize='small',fontproperties=font ,loc='right')



##################  False Alarms ######################

ax2=plt.subplot(2,2,3, projection=ccrs.PlateCarree())
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

dbz_levels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
dbz_rgb = np.array([[4,233,231],
                    [1,159,244], [3,0,244],
                    [2,253,2], [1,197,1],
                    [0,142,0], [253,248,2],
                    [229,188,0], [253,149,0],
                    [253,0,0], [212,0,0],
                    [188,0,0],[248,0,253],
                    [152,84,198]], np.float32) / 255.0
dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb,\
                                           extend="max")
contour_c=ax2.pcolormesh(false_lons, false_lats, false_rain,cmap=dbz_map,norm=dbz_norm)
#contour_c=plt.pcolormesh(lons,lats, accu_rain1,cmap=cmap_rain,
#                        norm=norm_rain)

#cb = plt.colorbar(contour_c, orientation="vertical",
#                      pad=0.02, aspect=16, shrink=0.8)
#cb.set_label('mm/hr',size=10)
#cb.ax.tick_params(labelsize=10)
#    plt.title("Precipitation Rate @ " + str(FHR) +'z'+  '_'+str(FHR_date))
font = FontProperties()
font.set_family('serif')
font.set_name('Times')
ax2.text(0.05,0.95, 'Total Cases:'+ str(len(falsecount_file)),transform=ax2.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
ax2.set_title("700mb Median Refl for False Alarms",fontsize=14,fontproperties=font )
#ax2.title('Core Files',fontsize='small',fontproperties=font ,loc='right')

########################### Correct Negatives ##########

#plt.subplot(2,2,2)
ax3=plt.subplot(2,2,2, projection=ccrs.PlateCarree())
ax3.coastlines(resolution="50m",linewidths=1)
ax3.set_extent([cn_lonc-3,cn_lonc+3,cn_latc-3,cn_latc+3])
#        print(lon_e[m])
#        print(lat_e[m])
gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([cn_lonc-3,cn_lonc,cn_lonc+3])
gl.ylocator = mticker.FixedLocator([cn_latc-3,cn_latc,cn_latc+3])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

dbz_levels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
dbz_rgb = np.array([[4,233,231],
                    [1,159,244], [3,0,244],
                    [2,253,2], [1,197,1],
                    [0,142,0], [253,248,2],
                    [229,188,0], [253,149,0],
                    [253,0,0], [212,0,0],
                    [188,0,0],[248,0,253],
                    [152,84,198]], np.float32) / 255.0
dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb,\
                                           extend="max")
contour_c=ax3.pcolormesh(cn_lons, cn_lats, cn_rain,cmap=dbz_map,norm=dbz_norm)
#contour_c=plt.pcolormesh(lons,lats, accu_rain1,cmap=cmap_rain,
#                        norm=norm_rain)

#cb = plt.colorbar(contour_c, orientation="vertical",
#                      pad=0.02, aspect=16, shrink=0.8)
#cb.set_label('mm/hr',size=10)
#cb.ax.tick_params(labelsize=10)
font = FontProperties()
font.set_family('serif')
font.set_name('Times')
ax3.text(0.05,0.95, 'Total Cases:'+str(len(cncount_file)),transform=ax3.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
ax3.set_title("700mb Median Refl for Correct Negatives",fontsize=14,fontproperties=font )
#qx3.title('Core Files',fontsize='small',fontproperties=font ,loc='right')

############################### Misses ################################

#plt.subplot(2,2,4)
ax4=plt.subplot(2,2,4, projection=ccrs.PlateCarree())
ax4.coastlines(resolution="50m",linewidths=1)
ax4.set_extent([miss_lonc-3,miss_lonc+3,miss_latc-3,miss_latc+3])
#        print(lon_e[m])
#        print(lat_e[m])
gl = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([miss_lonc-3,miss_lonc,miss_lonc+3])
gl.ylocator = mticker.FixedLocator([miss_latc-3,miss_latc,miss_latc+3])
gl.xformatter= LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size':10, 'color':'black'}
gl.ylabel_style = {'size':10, 'color':'black'}

dbz_levels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
dbz_rgb = np.array([[4,233,231],
                    [1,159,244], [3,0,244],
                    [2,253,2], [1,197,1],
                    [0,142,0], [253,248,2],
                    [229,188,0], [253,149,0],
                    [253,0,0], [212,0,0],
                    [188,0,0],[248,0,253],
                    [152,84,198]], np.float32) / 255.0
dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb,\
                                           extend="max")

contour_c=ax4.pcolormesh(miss_lons, miss_lats, miss_rain, cmap=dbz_map,norm=dbz_norm)
cb_dbz = fig.colorbar(contour_c, ax=[ax1,ax2,ax3,ax4], orientation="vertical",
                      pad=0.02, aspect=16, shrink=0.8)
cb_dbz.set_label('dbz',size=10)
cb_dbz.ax.tick_params(labelsize=10)
font = FontProperties()
font.set_family('serif')
font.set_name('Times')
ax4.text(0.05,0.95, 'Total Cases:'+str(len(misscount_file)),transform=ax4.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
ax4.set_title("700mb Median Refl for Misses",fontsize=14,fontproperties=font )

fig.savefig(picture_dir+'700mb_Median_Refl_multipanel_24hr'+'.png')
