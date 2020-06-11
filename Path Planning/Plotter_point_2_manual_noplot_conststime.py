import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style as style
style.use('seaborn-talk') #sets the size of the charts
style.use('ggplot')

Vdr = 15  #drone cruise speed
thov = 7  #hovering time 5s + 2*1s
ttot = 30*60  #total flight time
layers = 3  #total number of layers
redo = 10   #total redo of each point
tib =     #time (seconds) between measurement points


f = open("file.txt",'w')


def points(r,s):
    return int((ttot-2*r/Vdr + tib)/(thov+tib))

#
# for sp in np.linspace(0,2000,10):
#     for ra in np.linspace(0,2.5,10):
#         f.writelines(str(sp)+','+str(ra)+','+str(points(ra,sp))+'\n')
#f.close()


radm = 5#km
space = 200 #----------------------


xl = np.linspace(-radm,radm,int(radm*2*1000/space))
yl = np.linspace(-radm,radm,int(radm*2*1000/space))
allc = []
for xxl in xl:
    for yyl in yl:
        if yyl**2 + xxl**2 <= radm**2:
            allc.append([xxl,yyl,(yyl**2 + xxl**2)**0.5])


allc.sort(key=lambda x: x[2])
allc = allc[::-1]




t = []
d = []
print('Num of nodes ',len(allc))
while len(allc) > 0:
    radd = 1000*allc[0][2]
    pnts = points(radd,space)
    if pnts > len(allc):
        pnts = len(allc)
    t.append(pnts*thov) #hover
    t.append((pnts-1)*tib) #between
    d.append((pnts - 1) * space)
    t.append(2*radd/Vdr) #transfer
    d.append(radd*2)

    allc = allc[pnts:]
    print('Transfer distance: ',round(radd,2),' [m], points covered: ',round(pnts,2),', time taken: ',round(pnts*thov+(pnts-1)*tib+2*radd/Vdr,2),' [s], points left: ',len(allc),' distance covered: ',round(sum(d)/1000,2), 'km')
    d = []


print('\nTime [hrs] ',round(redo*layers*sum(t)/3600,2))


