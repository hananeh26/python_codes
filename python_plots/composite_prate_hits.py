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
from scipy import stats

names=[]
genesis=[]
not_genesis=[]
Hit= []
False_alarms=[]
Correct_negatives=[]
genesis_time=[]
not_genesis_time=[]
em=[]
not_genesis_FT=[]
Hit_genesis=[]
Hit_genesis_time=[]

dirs='/work/noaa/aoml-hafs1/galaka/noscrub/B220/'


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
################################ Find when wind speed when it hit 34kt in not genesis list #############################

for gfile in not_genesis:
    Gfile=pd.read_csv(dirs+gfile,delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=headers,index_col=False)
# read the files in not_genesis list

    for t in range(len(Gfile.Vmax)):
        if Gfile.Vmax[t]>=34:            # find out the first hour when wind speed is above 34 kt
            not_genesis_FT.append(Gfile.Forecast_Period[t])
            break

############################## Find hit, misses, false alarms, and correct negatives ##################################

Hit=0
False_alarms=0
regex = r'invest+[0-40]'
genesis.sort()
for k in range(len(genesis)):
    if re.search(regex,genesis[k]):
        Hit_genesis.append(genesis[k])
        Hit_genesis_time.append(genesis_time[k])
        Hit += 1
    else:
        False_alarms +=1
            # invest90..
Misses=0
Correct_negatives=0
for l in range(len(not_genesis)):
    if re.search(regex, not_genesis[l]):
        Misses += 1
    else:
        Correct_negatives +=1


for t in range(len(Hit_genesis)):
    FHR_date=int(Hit_genesis[t].split('.')[1][:10])
    my_frame = pd.DataFrame(data={'Files':[Hit_genesis[t]],'Model_initialization_time':[FHR_date],
                                   'Genesis_time_0hr':[Hit_genesis_time[t]],'Genesis_time_24hr':[int(Hit_genesis_time[t])+24]})

#    print(my_frame)


#invest27l.2020101712.trak.hwrf.atcfunix (hit_genesis)
check_hit=[]
check_hit_time=[]
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
        check_hit.append(Hit_genesis[w])
        check_hit_time.append(Hit_genesis_time[w])

#print(len(check_hit_time))
#print(len(check_hit))

 
 
####################################################################
sdir=u'/work/noaa/aoml-hafs1/hjafary/Dorian/dorian05l.2019082306/'
dires2=u'/work/noaa/aoml-hafs1/galaka/FOR_HANA/bdeck_invests/'
dires3=u'/work/noaa/aoml-hafs1/galaka/FOR_HANA/adeck_invests/'
picture_dir=u'/work/noaa/aoml-hafs1/hjafary/total_prec/'
####################################################################

grib_files='/work/noaa/aoml-hafs1/hjafary/hwrf_core/'
hwrf_atcf='/work/noaa/aoml-hafs1/galaka/noscrub/B220/'

slist=os.listdir(sdir)
swrf_out=[]
#accu_rain=[]
time=[]
invest_file=[]
bt_id=[]
hwrf_core=[]
adeck_list=[]

################## Grab all the directories with hwrfprs-core #########################

for (root, dirs, files) in os.walk(grib_files):
    for directory in dirs:
        if directory.endswith('L'):
            swrf_out.append(grib_files+directory+'/')

################# Read all the grib files in  direcotry ############################
swrf_out.sort()
#print(swrf_out)
#/work/noaa/aoml-hafs1/hjafary/hwrf_core/06L/
for lstfiles in swrf_out:
    for (root, dirs, files) in os.walk(lstfiles):
        for filename in files:
            if filename.endswith('.grb2'): #and int(filename.split('.')[0][-3:-1])<=40:
                hwrf_core.append(filename)

##########################################################################################

#adeck_list.sort()
new_list=[]
lon_e,lat_e=[],[]
hwrf_core.sort()
Hit_genesis.sort()
all_sid=[]
all_date=[]
all_fhr=[]
for at in range(len(check_hit)):
    all_date.append(str(check_hit[at].split('.')[1][:10])) #2019082306
    all_fhr.append(str(check_hit_time[at])) #48
    all_sid.append(str(check_hit[at].split('.')[0][-3:-1]))
#    storm_name=str(Hit_genesis[at].split('.')[0][:-3])

###################### HWRF_B ATCF ###########################
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

for x in range(len(FOUND_FILE)):
    my_frame = pd.DataFrame(data={'Files':[check_hit[x]],'Genesis_time_0hr':[check_hit_time[x]],'FOUND_FILE':[FOUND_FILE[x]]})


check_hit.sort()
for c in range(len(check_hit)):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    a_file=pd.read_csv(hwrf_atcf+check_hit[c],delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],
                            names=headers,index_col=False)   
    for j in range(0,len(a_file.Forecast_Period)-1):
#        if a_file.Forecast_Period[j+1] == a_file.Forecast_Period[j]:
#            continue
        if a_file.Forecast_Period[j+1] - a_file.Forecast_Period[j]==3 and \
               a_file.Forecast_Period[j] == int(str(check_hit_time[c]).zfill(3)): 
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
accu_rain1=[]
count_file=[]
FOUND_FILE.sort()

#for plot_list in new_list:
for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':{'shortName': 'prate','stepRange':str(FHR1)},'indexpath':''})

    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])


#    len_check=[]    
    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= lat_e[k]-2 and item  <= lat_e[k]+2:
            min_lat.append(index)
    lat_len.append(len(min_lat))


    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-2 and itm <= (lon_e[k])+2:
            min_lon.append(idx)
    lon_len.append(len(min_lon))

for k in range(len(FOUND_FILE)):
    FHR1 = int(FOUND_FILE[k].split('.')[-2][1:4])
    FHR_date1=int(FOUND_FILE[k].split('.')[1][:10])
    grbs = xr.open_dataset(grib_files+all_sid[k]+'L/'+all_date[k]+'/'+FOUND_FILE[k], engine='cfgrib',backend_kwargs={'filter_by_keys':\
               {'shortName': 'prate','stepRange':str(FHR1)},'indexpath':''})
   
    lats1=(grbs.variables['latitude'][:])
    lons1=(grbs.variables['longitude'][:])
 
    min_lat=[]
    for index,item in list(enumerate(lats1)):
        if item >= lat_e[k]-2 and item  <= lat_e[k]+2:
            min_lat.append(index)
    min_lon=[]
    for idx,itm in list(enumerate(lons1)):
        if itm >= (lon_e[k])-2 and itm <= (lon_e[k])+2:
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

#    print(b)
#    print(minlon)

    lons.append(lons1[minlon[:]])
    lats.append(lats1[minlat[:]])
    accu_rain1.append((grbs.variables['prate'][minlat[:],minlon[:]])*3600)
#    accu_rain1.append(rain)

    count_file.append(FOUND_FILE[k])
    lat_c.append(lat_e[k])
    lon_c.append(lon_e[k])


count_file.sort()
print(np.shape(accu_rain1))

z=np.nanmean(accu_rain1,axis=0)
print(np.shape(z))

#result=0
#for r in accu_rain1:
#    result=result+r
#z=(result/len(count_file))

lon_c=np.mean(lon_c)-360
lat_c=np.mean(lat_c)


#ly=0
#for item in lats:
#    ly=ly+item
#y=(ly/len(count_file))
#print(np.shape(y))
y=np.nanmean(lats,axis=0)

#lx=0
#for item in lons:
#    lx=lx+item
#x=(lx/len(count_file))
#print(np.shape(x))
x=np.nanmean(lons,axis=0)

fig = plt.figure(figsize=(10,7))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution="50m",linewidths=1)
ax.set_extent([lon_c-2,lon_c+2,lat_c-2,lat_c+2])
#        print(lon_e[m])
#        print(lat_e[m])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidths=1, color='black', linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([lon_c-2,lon_c,lon_c+2])
gl.ylocator = mticker.FixedLocator([lat_c-2,lat_c,lat_c+2])
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
contour_c=plt.pcolormesh(x, y, z, cmap=cmap_rain,norm=norm_rain)
#contour_c=plt.pcolormesh(lons,lats, accu_rain1,cmap=cmap_rain,
#                        norm=norm_rain)

cb = plt.colorbar(contour_c, orientation="vertical",
                      pad=0.02, aspect=16, shrink=0.8)
cb.set_label('mm/hr',size=10)
cb.ax.tick_params(labelsize=10)
#    plt.title("Precipitation Rate @ " + str(FHR) +'z'+  '_'+str(FHR_date))
font = FontProperties()
font.set_family('serif')
font.set_name('Times')
plt.title("Avg Precipitation Rate at Genesis",fontsize='small',fontproperties=font )
plt.title('Total Cases:'+ str(len(count_file)),fontsize='small',fontproperties=font ,loc='left')
plt.title('Core Files',fontsize='small',fontproperties=font ,loc='right')
plt.savefig(picture_dir+"Avg_Precipitation_Rate"+'.png')
plt.close()

