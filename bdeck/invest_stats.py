import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import re
import datetime
import glob,os

names=[]
genesis=[]
not_genesis=[]
#Hit= []
#False_alarms=[]
#Correct_negatives=[]
genesis_time=[]
not_genesis_time=[]
em=[]
not_genesis_FT=[]
Hit_genesis=[]
Hit_genesis_time=[]

dirs='/work/noaa/aoml-hafs1/galaka/noscrub/B220/'
names=[]
regex = r'invest+[0-99]+l'
for x in glob.glob('/work/noaa/aoml-hafs1/galaka/noscrub/B220/*.atcfunix'):
    if re.search(regex, x):
        names.append(os.path.basename(x))
    
#check=sorted([])
#for i in range(len(names)-1):
#    if names[i][10:14]!='2020':
#        check.append(dirs+names[i])

check=[]
for i in range(len(names)):
    check.append(dirs+names[i])
#print(check)

names.sort()
for sfile in sorted(names):
    headers = ['Basin','Invest_Number','Date_Time','TECHNUM','TECH','Forecast_Period','LatN','LonW','Vmax','MSLP','TY']
    A_file=pd.read_csv(dirs+sfile,delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=headers,index_col=False)
    
    FP=[]
    for i in range(len(A_file.Vmax)):
        if A_file.Vmax[i]>=34:
            FP.append(A_file.Forecast_Period[i])
#            print(A_file.Vmax)

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
                
################################# Categorize files into genesis and not genesis also append genesis time ##################


    if len(em)>=5:    # If the total len of em is greater than 5 --> [0,3,6,9,12], it meets genesis requirement. 
        genesis.append((sfile))
        genesis_time.append(em[0])
#        print(sfile,", 'Generated into TC' ")
#        print(genesis)
    else:
        not_genesis.append(sfile)

    em.clear()
        
#print(genesis)        
################################ Find when wind speed hit 34 in not genesis list ############################################
        
for gfile in sorted(not_genesis):
    Gfile=pd.read_csv(dirs+gfile,delimiter=',',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=headers,index_col=False)
# read the files in not_genesis list

    for t in range(len(Gfile.Vmax)):
        if Gfile.Vmax[t]>=34:            # find out the first hour when wind speed is above 34 kt
            not_genesis_FT.append(Gfile.Forecast_Period[t])
            break
    
############################## Find hit, misses, false alarms, and correct negatives #########################################
   
Hit=0
False_alarms=0
regex = r'invest+[0-40]'
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

print('Hits',Hit)
#print(len(Hit_genesis_time)) 
print('Correct_negatives',Correct_negatives)
print('False_alarms',False_alarms)
print('Misses',Misses)
print('sum',Misses+False_alarms+Correct_negatives+len(Hit_genesis))
############################### Print Statements #############################################################################

# print('Hit',Hit)
# print('False_alarms',False_alarms)
# print('Misses',Misses)
# print('Correct_negatives',Correct_negatives)
#    my_frame = pd.DataFrame(data={'Hit':[Hit],'False_alarms':[False_alarms],
#    'Misses':[Misses],'Correct_negatives':[Correct_negatives],
#    'Total':[Hit+False_alarms+Misses+Correct_negatives]})

for t in range(len(Hit_genesis)):
    FHR_date=int(Hit_genesis[t].split('.')[1][:10])
    my_frame = pd.DataFrame(data={'Files':[Hit_genesis[t]],'Model_initialization_time':[FHR_date],
                                   'Genesis_time_0hr':[Hit_genesis_time[t]],'Genesis_time_24hr':[int(Hit_genesis_time[t])+24]})

#    print(my_frame)
#    fig = plt.figure(figsize = (8, 2))
#    ax = fig.add_subplot(111)

#    ax.table(cellText = my_frame.values,
#                 colLabels = my_frame.columns,
#                 loc = "center")
#    ax.set_title("2017-2020 HWRF-B Forecast")

#    ax.axis("off");
#    plt.savefig('/work/noaa/aoml-hafs1/hjafary/bdeck/HWRF-B-Forecast',dpi=500)
#    plt.close()   


