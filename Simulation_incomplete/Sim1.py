# -*- coding: utf-8 -*-
'''

	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

   Created on Tue Jun  9 22:56:31 2020
	@author: moroberson, vgavra
   @version 2.0: fixed pygame.quit(), cleaned layout, modified plume model


	Project Supervisors:
		- Dr. Irene C. Dedoussi
		- Dr. Ir. Mirjam Snellen
		- Ir. Lorenzo Pasqualetto Cassinis
		- Mark Schelbergen

	This script is used for the initial sizing of the propellers / aerodynamic surfaces for the calcConcentrationept selection
'''
import numpy as np
import pygame as pg
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import math


# Constants:
    #colours
black     = (0,0,0)
white     = (255,255,255)
lightblue = (64,224,238)
blue      = (64,22,255)
red       = (255,0,0)
difblue   = (0,76,153)
lightred  = (255, 102, 102)

   # pygame screen
xmax = 1200
ymax = 800
sf   = 1 #scaling factor


# %%+++++++++++++++++++++++++++++++ Plume Modelling ++++++++++++++++++++++++++++++++++++++++++++++

class Plume():

    """ Class to generate the aircraft emission plume based on the Gaussian model"""

    def __init__(self,source=[10,10],thresh=0.001):

        """
        Constructor for the plume class
        """
        self.spnum  = 100            #resolution of plume - higher is better must more computationally expensive
        self.ym     = 30000*10          #  m-  width of plume taken into account from aircraft longitudinal axis
        self.xm     = 100000          # m --length of plume taken into account from plane
        self.Q      = 270.480176      # mg/m/3rough
        self.v      = 8               # m/s - rough
        self.thresh =   0.00000770      # 1/% of max concentration

    def calcConcentration(self,lc):
        """
         Calcualtes the concentration given the aircraft location coordinates
             Input: lc - tuple of aircraft coordinates (source of the plume)
             Output: concentration filed as sclar function of (x,y)
        """
        x = (lc[0])
        y = (lc[1])
        #sigy = (2*ky*x/v)**0.5
        #sigz = (2*kz*x/v)**0.5

        self.sigz = 0.06 *x *(1 + 0.0015*x)**(-0.5)
        self.sigy = max(self.sigz, 0.08*x*(1 + 0.0001*x)**(-0.5))

        return (self.Q/( 2* np.pi*self.sigy*self.sigz*self.v))*np.exp(-y**2/(2*self.sigy**2))  #dont forget the 2!


    # def noemConcentration(self, concentration):
    #     """
    #     Normalises the concentration value a.k.a. transforms the Plume(0,sigma) to N(0,1)
    #     """

    #     return cocentration / (self.Q/(np.pi*sigy*sigz*self.v))

    def plume(self, ac,rot): #rot in rad

        """

        """
        rot = -(np.pi-rot)#corrected for right x-positive as 0deg
        xr = np.linspace(0, self.xm, self.spnum)
        yr1 = np.linspace(0, -self.ym, self.spnum)
        yr2 = np.linspace(0, self.ym, self.spnum)


        moveon = False
        xlst1 = []
        xlst2 = []
        xlst = []
        ylst1 = []
        ylst2 = []
        ylst = []



        for xxl in xr:
            moveon = False
            for yyl in yr1:
                if self.calcConcentration((xxl, yyl)) < self.thresh and moveon == False:
                    if len(ylst1) > 0:
                        if yyl == ylst1[-1]:
                            moveon = True
                    # print('appended 1',xxl,yyl)
                    ylst1.append(yyl)
                    xlst1.append(xxl)
                    moveon = True

        for xxl in xr:
            moveon = False
            for yyl in yr2:
                if self.calcConcentration((xxl, yyl)) < self.thresh and moveon == False:
                    if len(ylst2) > 0:
                        if yyl == ylst2[-1]:
                            moveon = True
                    ylst2.append(yyl)
                    xlst2.append(xxl)
                    # print('appended 2',xxl,yyl)
                    moveon = True

        xlst2old = xlst2[:]

        #removing not needed tails
        for l in range(len(xlst1))[::-1]:
            if ylst1[l] == 0:
                xlst1.pop(l)
                ylst1.pop(l)
        for l in range(len(xlst2))[::-1]:
            if ylst2[l] == 0:
                xlst2.pop(l)
                ylst2.pop(l)

        #combining and adding the y=0 part
        xlst = xlst2 + [xlst2old[len(xlst2)]] + xlst1[::-1]
        ylst = ylst2 + [0] + ylst1[::-1]


        flst = []
        for i in range(len(xlst)):
            xlr = (xlst[i] * np.cos(rot) - ylst[i] * np.sin(rot)) + ac[0] #tempplume
            ylr = (xlst[i] * np.sin(rot) + ylst[i] * np.cos(rot)) + ac[1] #tempplume
            flst.append((xlr,ylr))

        return flst


# %% +++++++++++++++++++++++++++++++ Mission Geometry ++++++++++++++++++++++++++++++++++++++++++++++

layers = 4
laymin = 10 #m
laystep = 10 #m
radm = 5000  # m  (max radius)
space = 150  #m

vps = 80/1000 #ac speed m/(0.001*s)
vpsfix = 0.99 #factor on vps
vps = vps*vpsfix



#LTO
lastLT0takeoff = - 40 - 100
lastLTOlanding = - 30 - 22
tLTO= 50
LTOf = 0.1 #LTO factor reduce tLTO

tf = 1 #timefactor 1=realtime
#funcs


# %% +++++++++++++++++++++++++++ Pygame Helper Functions ++++++++++++++++++++++++++++++++++++++++++
#distance m to pixel
def d(dist): #length to pixel
    return int(dist*(1/16))

def tdc(c): #coords to pixel
    return (int(d(c[0])+xmax/2),int(-d(c[1])+ymax/2))

def tdcl(cl,xscreen=xmax,yscreen=ymax):
    """
    Returns the converted list of coordinates to pixels
    Input:  cl = list of coordinates
            xcreen = x-dimension of the screen
            yscreen =  y-dimension of the screen
    """

    lstc = []
    for c in cl:
        lstc.append((int(d(c[0])+xscreen/2),int(-d(c[1])+yscreen/2)))

    return lstc

def getPos():
    pos = pg.mouse.get_pos()
    print(pos)
    return (pos)

def lin(pos,x1,y1,x2,y2):
    return ((y2-y1)/(x2-x1))*(pos-x1)+y1

def vdir(d,P1,P2):
    x1,y1 = P1
    x2,y2 = P2
    return (x1+(x2-x1)*d,y1+(y2-y1)*d)




#%% +++++++++++++++++++++++++++++++++ Airport Data +++++++++++++++++++++++++++++++++++++++++++++++++
#define no go zones:
rwy36Lcoords = [(-3482,2307), (-2882,12060), (-2670,12000), (-3270,2247)]
rwy24coords = [(1424,-245),(1600,-500),(-4500,-3890),(-5000,-6000),(-6000,-4409)]
rwy36Lpoly = Polygon(rwy36Lcoords)
rwy24poly = Polygon(rwy24coords)

rwy36Lstart = np.array(((rwy36Lcoords[0][0]+rwy36Lcoords[3][0])/2,(rwy36Lcoords[0][1]+rwy36Lcoords[3][1])/2))
rwy36Lfinish = np.array(((rwy36Lcoords[1][0]+rwy36Lcoords[2][0])/2,(rwy36Lcoords[1][1]+rwy36Lcoords[2][1])/2))
rwy24start = np.array(((rwy24coords[0][0]+rwy24coords[1][0])/2,(rwy24coords[0][1]+rwy24coords[1][1])/2))
dxc = (rwy24coords[0][0]-rwy24start[0])
dyc = (rwy24coords[0][1]-rwy24start[1])
rwy24mid = np.array(((rwy24coords[2][0]+dxc),(rwy24coords[2][1]+dyc)))
rwy24finish1 = np.array(((rwy24coords[4][0]-dxc),(rwy24coords[4][1]-dyc)))
rwy24finish2 = np.array(((rwy24coords[3][0]+dxc),(rwy24coords[3][1]+dyc)))

rwy36Langle = math.atan2(rwy36Lfinish[1]-rwy36Lstart[1],rwy36Lfinish[0]-rwy36Lstart[0])
rwy24angle = math.atan2(rwy24start[1]-rwy24finish1[1],rwy24start[0]-rwy24finish1[0])


#plume (tempplume)
# tempplumecoords = tdcl(Plume.plume((0,0),rwy24angle))


#nodes
xl = np.linspace(-radm, radm, int(radm * 2 / space)+1)
yl = np.linspace(-radm, radm, int(radm * 2 / space)+1)
zl = np.linspace(laymin, laystep * layers + laymin, layers)
coords = []
for xxl in xl:
    for yyl in yl:
        if rwy36Lpoly.contains(Point(xxl,yyl)) == False and rwy24poly.contains(Point(xxl,yyl)) == False:
            if yyl ** 2 + xxl ** 2 <= radm ** 2:
                coords.append([xxl, yyl, (yyl ** 2 + xxl ** 2) ** 0.5])



#aircraft

actypes = np.array(['A300','A318','A319','A320','A321','A330','A340','A350','A380','AT72','B737',
                     'B747','B757','B767','B777','B787','A220','E170','E190','CRJ700'])
acfreq = np.array([11,6,33,77,22,33,8,16,13,10,89,25,7,25,31,44,12,11,14,7])#initial guess
acfreq = acfreq/(sum(acfreq))#normalised sum -> 1
aclivetakeoff = []
acliveland = []
actakeoffpos = []
aclandpos = []
actot = []
accol = []
for px in range(len(actypes)):

    mx = len(actypes)/255
   #print((px/mx,255-px/mx,abs(256/2 - px/mx)))
    accol.append((px/mx,255-px/mx,abs(256/2 - px/mx)))




# %% ++++++++++++++++++++++++ Main Simulation +++++++++++++++++++++++++++++++++++++++++++++++++++++
#Start of Pygame
scr = pg.display.set_mode((xmax*sf,ymax*sf))
win = pg.Surface((xmax, ymax)) #leave

pg.init()

#Sim Loop
t = pg.time.get_ticks()*0.001

#initials

running = True
while running == True:
    t = tf*t
    pg.event.pump()
    t0 = t
    t = pg.time.get_ticks()*0.001
    dt = t - t0




    scr.fill(white)

    #5km radius
    pg.draw.circle(scr,red,(int(xmax/2),int(ymax/2)),d(5000),d(50))

    #geocaged
    pg.draw.polygon(scr,lightred,tdcl(rwy36Lcoords),d(40))
    pg.draw.polygon(scr, lightred,tdcl(rwy24coords),d(40))


    #coords of nodes
    for coord in coords:
        pg.draw.circle(scr,lightblue,tdc(coord),2)



    #centre node
    pg.draw.circle(scr, red, tdc((0,0)), 5)

    #flight path points
    for l in [rwy36Lstart,rwy36Lfinish,rwy24start,rwy24finish1,rwy24finish2,rwy24mid]:
        pg.draw.circle(scr, blue, tdc(l), 2)


    #aircraft

    #T/O
    if lastLT0takeoff + 2*tLTO*LTOf < t:          #add some noise (random........ ) !!!!!!!!!!!!!!!!!
        act = np.random.choice(actypes,p=acfreq)
        accolval = accol[(list(actypes)).index(act)]
        aclivetakeoff.append([act,rwy36Lstart,t,rwy36Langle,accolval])  # [type,pos,time added,angle,colour]
        print('Take-off -- added: ',act)#make random ############ by freq
        actot.append(act)
        lastLT0takeoff = t

    for ac in aclivetakeoff:
        pg.draw.circle(scr, ac[4], tdc(ac[1]) , 3)
        ac[1] = vdir((t-ac[2]+0.001)*vps,rwy36Lstart,rwy36Lfinish)

        if rwy24poly.contains(Point(ac[1])) == False and rwy36Lpoly.contains(Point(ac[1])) == False:
            aclivetakeoff.remove(ac)
            print('Take-off -- removed: ',str(ac[0]))




    # Landing
    if lastLTOlanding + 2*tLTO*LTOf < t:  # add some noise (random........ ) !!!!!!!!!!!!!!!!!
        act = np.random.choice(actypes,p=acfreq)
        accolval = accol[(list(actypes)).index(act)]
        acliveland.append([act,rwy24start,t,rwy24angle,accolval])
        print('Landing -- added: ',act)
        actot.append(act)
        lastLTOlanding = t

    for ac in acliveland:
        pg.draw.circle(scr, ac[4], tdc(ac[1]) , 3)
        ac[1] = vdir((t-ac[2]+0.001)*vps,rwy24finish1,rwy24start)

        if rwy24poly.contains(Point(ac[1])) == False and rwy36Lpoly.contains(Point(ac[1])) == False:
            acliveland.remove(ac)
            print('Landing -- removed: ',str(ac[0]))




    #draw plumes
    for ac in acliveland+aclivetakeoff:
        acPlume = Plume()
        pg.draw.polygon(scr, lightred, tdcl(acPlume.plume(ac[1],ac[3])),d(40))



    pg.display.flip()



    #quit
    for event in pg.event.get():
        if event.type==pg.QUIT:
            running = False

    keys = pg.key.get_pressed()
    if keys[pg.K_ESCAPE]:
        running = False


pg.display.quit()
pg.quit()

print('\n\n')
actots = list(set(actot))
actots.sort()
for acs in actots:
    print(acs,' -- count: ',actot.count(acs),', freq: ',round(actot.count(acs)*100/len(actot),3),'%')

print(" \n \n Done! ")

