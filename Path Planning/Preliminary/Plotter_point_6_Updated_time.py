import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib.style as style

Vdr = 10
thov = 7
ttot = 30*60
layers = 4
redo=10

import matplotlib
import matplotlib.pyplot as plt


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







datapnts = 7 #number of datapoints
tmeas = 5#sec measuring

redo = 3 #redo eachpoint   leave out
#spacescale = 0.9
dsp = 1
Fset = 2/7 #per day
Fminn = 1/25###############################


nc = 0.9 #navigation contingency
rex = 0.6 #for any manouvering extras and SAA/DAA
redo = redo + rex

laymin = 10 #m
laystep = 10 #m
tmeaspday = 18 #operational hours per day  (noncurfew+1hr)

actrwy = 2 #number of active runways whether departing or landing
radm = 5  # km  (max radius)

LTOtbetween = 50 #s time between LTO
airportoptime = 20 #hrs
acLTO = 100 #LTO event time     # 87m/s for 100sec = 8.7km



def points(r,s):
    return int((ttot-2*r/Vdr + s/Vdr)/(thov+s/Vdr))

#f = open("file.txt", 'w')
# for sp in np.linspace(0,2000,10):
#     for ra in np.linspace(0,2.5,10):
#         f.writelines(str(sp)+','+str(ra)+','+str(points(ra,sp))+'\n')
#
# f.close()

radm = 5  # km
space = 150

rlst = []
pcovlst = []
ttlst = []


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

    print('radius: ',round(radd,2),', points covered: ',round(pnts,2),', time taken: ',round(pnts*thov+(pnts-1)*space/Vdr+2*radd/Vdr,2),', len left: ',len(allc))
    rlst.append(radd)
    pcovlst.append(pnts)
    ttlst.append(pnts*thov+(pnts-1)*space/Vdr+2*radd/Vdr)
print('\nSpacing: ',space,'  --  time [hrs] ',round(layers*sum(t)/3600,4))


plt.scatter(rlst,ttlst,s=11,color=(0,0.5,0.9))

plt.axvline(0, color=(0, 0, 0), linewidth=0.8)
#plt.axhline(0, color=(0, 0, 0), linewidth=0.8)
plt.xlabel('Distance from Ground Station [m]')
plt.ylabel('Flight Time [s]')
plt.minorticks_on()  # set minor ticks
plt.grid(which='major', linestyle='-', linewidth='0.3', color='black')  # customise major grid
plt.grid(which='minor', linestyle=':', linewidth='0.3', color='grey')  # customise minor grid
plt.axhline(ttot,linestyle='--',color='r')

plt.show()