import numpy as np
import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as style


texpsize = [15, 18, 22]

plt.rc('font', size=texpsize[1], family='serif')  # controls default text sizes
plt.rc('axes', titlesize=texpsize[1])  # fontsize of the axes title
plt.rc('axes', labelsize=texpsize[1])  # fontsize of the x and y labels
plt.rc('xtick', labelsize=texpsize[0])  # fontsize of the tick labels
plt.rc('ytick', labelsize=texpsize[0])  # fontsize of the tick labels
plt.rc('legend', fontsize=texpsize[0])  # legend fontsize
plt.rc('figure', titlesize=texpsize[2])  # fontsize of the figure title
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['axes.facecolor'] = 'white'
matplotlib.rcParams["legend.fancybox"] = False






#old
# style.use('seaborn-talk') #sets the size of the charts
# style.use('ggplot')


datapnts = 7 #number of datapoints
tmeas = 5#sec measuring

redo = 3 #redo eachpoint   leave out
spacescale = 0.994
Fset = 2/7 #per day
Fminn = 1/20###############################


nc = 0.9 #navigation contingency

Vdr = 10
thov = 7
ttot = 30*60
layers = 4
laymin = 10 #m
laystep = 20 #m
tmeaspday = 17 #operational hours per day

actrwy = 2 #number of active runways whether departing or landing
radm = 5  # km  (max radius)

LTOtbetween = 50 #s time between LTO
airportoptime = 20 #hrs
acLTO = 100 #LTO event time     # 87m/s for 100sec = 8.7km

actype = ['A','B','C']
acfreq = [4,3,55] #freq (per day)


if airportoptime < tmeaspday:
    print('Airport operational time is less than measuring time')
    sys.exit()

if airportoptime*3600 < sum(acfreq)*LTOtbetween/actrwy:
    print('Aircraft LTO not consistent with airport operation time')
    sys.exit()







def points(r,s):
    return int((ttot-2*r/Vdr + s/Vdr)/(thov+s/Vdr))





def Flim(space,opdays):

    optime = 24 * opdays * (tmeaspday / 24)  # operational time total hrs (adjusted for tmeaspday)
    optime = optime * 3600  # convert to seconds
    # FINDING POINTS WITHIN CYLINDER RADIUS ----
    xl = np.linspace(-radm, radm, int(radm * 2 * 1000 / space))
    yl = np.linspace(-radm, radm, int(radm * 2 * 1000 / space))
    zl = np.linspace(laymin,laystep*layers+laymin,layers)
    allc0 = []
    for xxl in xl:
        for yyl in yl:
            if yyl ** 2 + xxl ** 2 <= radm ** 2:
                allc0.append([xxl, yyl,(yyl ** 2 + xxl ** 2)**0.5])
    #Sort for radius
    allc0.sort(key=lambda x: x[2])  #2 for radius 3 for score
    allc0 = allc0[::-1]
    #add z_height per layer
    allc = []
    for zzl in zl:
        for i in allc0:
            allc.append([i[0],i[1],zzl,i[2]])
    #now sorted for each z layer and then decreasing radius within


    #RUNNING THROUGH NODES FOR TIME
    t = []
    tusetemp = [] #useful measuring time (temp)
    tuse = [] #useful measuring time
    twaste = []
    allctemp = allc #for 2nd while
    tuseday = [] #each day measured time
    day = 1
    while sum(t)*nc < optime:
        while len(allctemp) > 0 and sum(t)*nc < optime:
            radd = 1000*allctemp[0][3]
            pnts = points(radd,space)  #Assuming layer transition negligible
            if pnts > len(allctemp): #if at last points
                pnts = len(allctemp)
            t.append(pnts*thov) #hover
            tusetemp.append(pnts*thov) #hover (measuring)
            tuse.append(pnts*thov)
            t.append((pnts-1)*space/Vdr) #between
            twaste.append((pnts-1)*space/Vdr)
            t.append(2*radd/Vdr) #transfer
            twaste.append(2*radd/Vdr)
            allctemp = allctemp[pnts:]
            mfrac = sum(tuse) / sum(twaste)
            if sum(t)*nc > day*60*60*24*(tmeaspday/24):
                tuseday.append(sum(tusetemp))
                mfrac = sum(tuse) / (sum(twaste)+sum(tuse))
                #print('Day ', day, ': ', round(sum(tusetemp) / 3600, 3), ' hrs of measurements')
                tusetemp = []
                day += 1



        if len(allctemp) <= 0: #restarts while loop
            allctemp = allc


    #print('\nSpacing: ',round(space,2),'  mfrac: ',round(100*mfrac,2),'%')
    Flimval = datapnts*tmeas/(LTOtbetween*mfrac)
    return Flimval/(optime/3600/24) #returns limiting freq [flights/day]


def tm(space):

    #space = 200  # ----------------------
    xl = np.linspace(-radm, radm, int(radm * 2 * 1000 / space))
    yl = np.linspace(-radm, radm, int(radm * 2 * 1000 / space))
    allc = []
    for xxl in xl:
        for yyl in yl:
            if yyl ** 2 + xxl ** 2 <= radm ** 2:
                allc.append([xxl, yyl, (yyl ** 2 + xxl ** 2) ** 0.5])


    allc.sort(key=lambda x: x[2])
    allc = allc[::-1]

    t = []
    #print(len(allc))
    while len(allc) > 0:
        radd = 1000*allc[0][2]
        pnts = points(radd,space)
        if pnts > len(allc):
            pnts = len(allc)
        t.append(pnts*thov) #hover
        t.append((pnts-1)*space/Vdr) #between
        t.append(2*radd/Vdr) #transfer
        allc = allc[pnts:]

        #print('radius: ',round(radd,2),', points covered: ',round(pnts,2),', time taken: ',round(pnts*thov+(pnts-1)*space/Vdr+2*radd/Vdr,2),', len left: ',len(allc))

    #print('\nSpacing: ',space,'  --  time [hrs] ',round(layers*sum(t)/3600,4))
    return round(redo*layers*sum(t)*nc/3600/24,3)



space = 300
opdays = 30 #days


Flimval = 2

spacelst = []
timelst = []
Flimlst = []

while Flimval > Fminn:
    opdays = tm(space)
    Flimval = Flim(space, opdays)

    opfact =  (tmeaspday / 24)
    print('\n----------------------------')
    print('spacing: ',space)
    print('time: ',opdays/opfact)#   real time: ',tm(space))
    print('F_lim: ',Flimval)
    print('difference: ',-(Fset - Flimval)/Fset)
    spacelst.append(space)
    timelst.append(opdays/opfact)
    Flimlst.append(7*Flimval)
    space = space * spacescale

    #print(space,Fset,Flimval,-(Fset - Flimval)/Fset)




c1 = (0,102/256,204/256)
c2 = (1,0,0)



plt.axvline(0,color=(0,0,0),linewidth=0.8)
plt.axhline(0,color=(0,0,0),linewidth=0.8)

plt.minorticks_on() # set minor ticks
plt.grid(which='major', linestyle='-', linewidth='0.3', color='black') # customise major grid
plt.grid(which='minor', linestyle=':', linewidth='0.3', color='grey') # customise minor grid
plt.axvline(7*Fset,color=(0,153/256,76/256),linestyle='--')


plt.plot(Flimlst,spacelst, color=c1)

plt.ylabel('Mesh Spacing [m]')
plt.xlabel('Limiting A/C Type Freq [1/week]')


plt.show()