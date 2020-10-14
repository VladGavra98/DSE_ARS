import numpy as np
import shapely as sh
from shapely.geometry import Point,LineString
from shapely.geometry.polygon import Polygon
import matplotlib.pyplot as plt
from collections import deque, namedtuple
import six
import sys
sys.modules['sklearn.externals.six'] = six
import mlrose
import random
import time

start_time = time.time()

print('\n\nInitialising\n\n------\n')
radm = 5000
space = 550
sfac = 0.17  #size of circles and linewidth
radfrac = 0.0 #fraction of maxradius below which points are not checked for dijkstra #0 is no effect

Vdr = 15
thov = 7
ttot = 36*60

#Ground station moved
gsx = -3500 #xlocation [m]
gsy = 500 #ylocation [m]


#Fixed -> for visual once routes known:
fix = False  #True means pre-routed as already determined when set to False
f550 = False #550m fix for visualisation purposes #default = False


#444m
#routesfix = [38, 172, 296, 37, 111, 276, 368, 11, 158, 51, 249]#444
#routesfixc = [137, 220, 263, 46, 167, 209, 269, 141, 265, 0, 203]#444

#550
#routesfix = [22, 114, 214, 247, 2, 10, 46, 198, 138]#550
#routesfixc = [102, 137, 172, 156, 69, 106, 125, 0, 138]#550
#f550 = True

#600
#routesfix = [20, 105, 19, 59, 148, 176, 28]#600  #furthest point routes
#routesfixc = [93, 139, 36, 84, 98, 113, 0]#600 #closest point routes


#------------------------------------------------------------------------------------------------#
#Dijkstra taken from https://dev.to/mxl/dijkstras-algorithm-in-python-algorithms-for-beginners-dkc
# we'll use infinity as a default distance to nodes.
inf = float('inf')
Edge = namedtuple('Edge', 'start, end, cost')

def make_edge(start, end, cost=1):
  return Edge(start, end, cost)


class Graph:
    def __init__(self, edges):
        # let's check that the data is right
        wrong_edges = [i for i in edges if len(i) not in [2, 3]]
        if wrong_edges:
            raise ValueError('Wrong edges data: {}'.format(wrong_edges))

        self.edges = [make_edge(*edge) for edge in edges]

    @property
    def vertices(self):
        return set(
            sum(
                ([edge.start, edge.end] for edge in self.edges), []
            )
        )

    def get_node_pairs(self, n1, n2, both_ends=True):
        if both_ends:
            node_pairs = [[n1, n2], [n2, n1]]
        else:
            node_pairs = [[n1, n2]]
        return node_pairs

    def remove_edge(self, n1, n2, both_ends=True):
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        edges = self.edges[:]
        for edge in edges:
            if [edge.start, edge.end] in node_pairs:
                self.edges.remove(edge)

    def add_edge(self, n1, n2, cost=1, both_ends=True):
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        for edge in self.edges:
            if [edge.start, edge.end] in node_pairs:
                return ValueError('Edge {} {} already exists'.format(n1, n2))

        self.edges.append(Edge(start=n1, end=n2, cost=cost))
        if both_ends:
            self.edges.append(Edge(start=n2, end=n1, cost=cost))

    @property
    def neighbours(self):
        neighbours = {vertex: set() for vertex in self.vertices}
        for edge in self.edges:
            neighbours[edge.start].add((edge.end, edge.cost))

        return neighbours

    def dijkstra(self, source, dest):
        assert source in self.vertices, 'Such source node doesn\'t exist'
        distances = {vertex: inf for vertex in self.vertices}
        previous_vertices = {
            vertex: None for vertex in self.vertices
        }
        distances[source] = 0
        vertices = self.vertices.copy()

        while vertices:
            current_vertex = min(
                vertices, key=lambda vertex: distances[vertex])
            vertices.remove(current_vertex)
            if distances[current_vertex] == inf:
                break
            for neighbour, cost in self.neighbours[current_vertex]:
                alternative_route = distances[current_vertex] + cost
                if alternative_route < distances[neighbour]:
                    distances[neighbour] = alternative_route
                    previous_vertices[neighbour] = current_vertex

        path, current_vertex = deque(), dest
        while previous_vertices[current_vertex] is not None:
            path.appendleft(current_vertex)
            current_vertex = previous_vertices[current_vertex]
        if path:
            path.appendleft(current_vertex)
        return path
#------------------------------------------------------------------------------------------------#



def intersect(a, b, p, q):
    '''Returns if intersects'''
    m2 = (p[1] - q[1]) / (p[0] - q[0])#never infty
    c2 = p[1] - m2 * p[0]
    if (a[0] - b[0]) == 0:#exception case gradient = infinity
        #print('made exception')
        xi = a[0]
        yi = m2*xi+c2
        if p[0] < xi < q[0] or q[0] < xi < p[0]:
            if a[1] < b[1]:
                if a[1] < yi < b[1]:
                    return True
            if a[1] > b[1]:
                if a[1] > yi > b[1]:
                    return True
    else:
        m1 = (a[1]-b[1])/(a[0]-b[0])
        c1 = a[1] - m1*a[0]

        xi = (c2-c1)/(m1-m2)

        if max(min(a[0],b[0]),min(p[0],q[0])) < xi < min(max(a[0],b[0]),max(p[0],q[0])):
            return True
    return False


def intersectpoly(p1,p2):
    '''Returns if a line (/route) between two points intersects geocaged runway'''
    hit = False
    for r36 in range(len(rwy36Lcoords)):
        if intersect(p1,p2,rwy36Lcoords[r36%len(rwy36Lcoords)],rwy36Lcoords[(r36+1)%len(rwy36Lcoords)]) == True:
            hit = True
    for r24 in range(len(rwy24coords)):
        if intersect(p1,p2,rwy24coords[r24%len(rwy24coords)],rwy24coords[(r24+1)%len(rwy24coords)]) == True:
            hit = True
    return hit



def uniq(lst):
    '''Part 1 uniquifies a list of lists'''
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item
def set2(l):
    '''Part 2'''
    return list(uniq(sorted(l, reverse=True)))

def points(r):
    '''Returns number of points that may be covered depending on distance'''
    return int((ttot-2*r/Vdr + space/Vdr)/(thov+space/Vdr))




def furthest2(coords,graphx): #furthest with constraints
    '''Returns furthest point (as route) and route with constrain on points looked at ( > radial value)'''
    mxval = []
    for i in coords:
        mxval.append((i[0]**2+i[1]**2)**0.5)
    maxrad = max(mxval)*radfrac-100
    ln = []
    #print(coords)
    lx = len(coords)
    print('Furthest Point Dijkstra -- 0%')
    for g,l in enumerate(coords):
        print('FPF: ',round(100*g/lx,1),'%') #Furthest Point Finding
        if (l[0]**2+l[1]**2)**0.5 > maxrad: #limits routes to be found that are wasted (assumes outermost point will ...
            rte = graphx.dijkstra(coordsorig[0][2], (l[2]))

            dist = []
            for i in range(len(rte)-1):

                dist.append(((coordsorig[rte[i]][0] - coordsorig[rte[i+1]][0])**2+(coordsorig[rte[i]][1] - coordsorig[rte[i+1]][1])**2)**0.5)
            ln.append(sum(dist))
    for i in range(len(ln)):
        if ln[i] == max(ln):
            print('Furthest Point Dijkstra -- 100%\n')
            return coords[i],max(ln)

def closest2(coords,graphx): #furthest with constraints
    '''Returns furthest point (as route) and route with constrain on points looked at ( > radial value)'''
    mxval = []
    coords2=[]
    for p in coords:
        coords2.append([coordsorig[p][0],coordsorig[p][1],p])
    coords = coords2[:]
    for i in coords:
        mxval.append((i[0]**2+i[1]**2)**0.5)
    maxrad = max(mxval)*radfrac-100
    ln = []
    #print(coords)
    lx = len(coords)
    print('Minimum Point Dijkstra -- 0%')
    for g,l in enumerate(coords):
        print('MPF: ',round(100*g/lx,1),'%') #Furthest Point Finding
        if (l[0]**2+l[1]**2)**0.5 > maxrad: #limits routes to be found that are wasted (assumes outermost point will ...
            rte = graphx.dijkstra(coordsorig[0][2], (l[2]))
            dist = []
            for i in range(len(rte)-1):

                dist.append(((coordsorig[rte[i]][0] - coordsorig[rte[i+1]][0])**2+(coordsorig[rte[i]][1] - coordsorig[rte[i+1]][1])**2)**0.5)
            ln.append(sum(dist))
    for i in range(len(ln)):
        if ln[i] == min(ln):
            print('Minimum Point Dijkstra -- 100%\n')
            return coords[i],min(ln)



def finddist(lst):
    '''Returns distance of a path (as list)'''
    dist = []
    rte = lst
    for i in range(len(rte) - 1):
        dist.append(((coords[rte[i]][0] - coords[rte[i + 1]][0]) ** 2 + (
                    coords[rte[i]][1] - coords[rte[i + 1]][1]) ** 2) ** 0.5)
    return sum(dist)



#Runway coordinates
rwy36Lcoords = [(-1000, -1700), (-1001 - 5000, 11700), (-1001 + 212 - 5000, 11701),
                (-1000 + 212, -1701)]  # no direct veritcal or horizotnal lines
rwy24coords = [(0, +1000), (212, 1000), (212 + 5000, -11701), (5001, -11700), (5000, -11701)]
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


#Nodes
xl = np.linspace(-radm, radm, int(radm * 2 / space)+1)
yl = np.linspace(-radm, radm, int(radm * 2 / space)+1)


coords = [[gsx,gsy,0]]
count = 1

for xxl in xl:
    for yyl in yl:
        if rwy36Lpoly.contains(Point(xxl,yyl)) == False and rwy24poly.contains(Point(xxl,yyl)) == False:
            if yyl ** 2 + xxl ** 2 <= radm ** 2:
                coords.append([xxl, yyl,count])
                count += 1




print('\nThere are ',len(coords),' nodes\n')


#Finding nodes that can be moved to for each node
cnext = []
cnextnum = [] #coord number (index) only
for loc in coords:

    tch = []
    tchnum = [] #coord number (index) only
    for check in coords:
        if 0.5 < ((loc[0]-check[0])**2 + (loc[1]-check[1])**2)**0.5 < 2**0.5*space*1.2:
            if intersectpoly(loc,check) == False:
                tch.append(check)
                tchnum.append(check[-1])
    cnext.append(tch)
    cnextnum.append(tchnum)



cx = 0

#Performing path loops to centre
furthestpoint = []
furthestroute = []
closestpoint = []
closestroute = []
run = True
removedpointsall = []
coordsorig = coords[:]
fcheck = [-2,-1]

if fix == True:
    furthestpoint = routesfix[:]
    print('---Using predetermined nodes---')


while len(coords) > 0 and fcheck[-1] != fcheck[-2] and run == True:
    coordsx = coords
    if fix == True:
        coordsx = [coordsorig[furthestpoint[cx]]]



    print('Making initial distance matrix (neighbour only)')
    lstalli = []
    for a,b,i in coords:
        for j in cnextnum[i]:
            ds = ((coordsorig[i][0] - coordsorig[j][0]) ** 2 + (coordsorig[i][1] - coordsorig[j][1]) ** 2) ** 0.5
            lstalli.append(((coordsorig[i][2]), (coordsorig[j][2]), ds))
    lstli = False  # if lstalli changed
    print('Made neighbouring nodes')
    graphi = Graph(lstalli)  # redo
    print('\nFinding Furthest Coord')
    try:
        coordfurthest,dsmax = furthest2(coordsx,graphi)
    except:
        print('!!Allowing old nodes to be traversed!!') #same coords (coords[]) but all routes possible (graphi)
        lstalli = []
        for a, b, i in coordsorig:
            for j in cnextnum[i]:
                ds = ((coordsorig[i][0] - coordsorig[j][0]) ** 2 + (coordsorig[i][1] - coordsorig[j][1]) ** 2) ** 0.5
                lstalli.append(((coordsorig[i][2]), (coordsorig[j][2]), ds))
        graphi = Graph(lstalli)  # redo
        coordfurthest, dsmax = furthest2(coordsx, graphi)
        lstli = True
    if dsmax == 0:
        print('!!V2 - Allowing old nodes to be traversed!! (d=0)')
        lstalli = []
        for a, b, i in coordsorig:
            for j in cnextnum[i]:
                ds = ((coordsorig[i][0] - coordsorig[j][0]) ** 2 + (coordsorig[i][1] - coordsorig[j][1]) ** 2) ** 0.5
                lstalli.append(((coordsorig[i][2]), (coordsorig[j][2]), ds))
        graphi = Graph(lstalli)  # redo
        coordfurthest, dsmax = furthest2(coordsx, graphi)

        lstli = True

    cx += 1

    print('Furthest Coord: ',coordfurthest,' Distance: ',dsmax,' [m]')

    coordfurthestroute = graphi.dijkstra(coordsorig[0][2], coordfurthest[2])
    print('Routing: ',coordfurthestroute)

    furthestroute.append(coordfurthestroute)
    if fix == False:
        furthestpoint.append(coordfurthest[2])
    pnts = points(dsmax)
    print('\nPoints to cover: ',pnts)
    if pnts < 1:
        print('Not enough points to measure!!')
    if pnts > len(coords):
        print('Need to adjust number of measurable points to points left')
        pnts = len(coords)

    lsttmp = []

    i = coordfurthest[2]

    print('\nCreating nearby distance matrix')
    nearbylst = [[0,0],[i,1]] #with station
    nearbylst = [[i, 0]]  # no station
    lstv = []
    for aa,bb,j in coords:
        if j != i and intersectpoly(coordsorig[i],coordsorig[j]) == False:
            ds = ((coordsorig[i][0] - coordsorig[j][0]) ** 2 + (coordsorig[i][1] - coordsorig[j][1]) ** 2) ** 0.5
            nearbylst.append([j,ds])


    nearbylst.sort(key=lambda x: x[1])
    nearbylst = nearbylst[:pnts]
    for i in nearbylst:
        lstv.append(i[0])

    # for ll,d in nearbylst:
    #     plt.scatter([coords[ll][0]],[coords[ll][1]])

    for l in coordsorig:
        if l[2] == None:
            print(l)


    lstall = []
    for j,x in nearbylst:
        for v in cnextnum[j]:
            if v in lstv:
                ds = ((coordsorig[v][0] - coordsorig[j][0]) ** 2 + (coordsorig[v][1] - coordsorig[j][1]) ** 2) ** 0.5
                if ds > 0:

                    if (coordsorig[j][2],coordsorig[v][2], ds) not in lstall:
                        lstall.append((coordsorig[v][2],coordsorig[j][2], ds))





    #plots all nearbys
    for l,s,r in lstall:  #would be nice
         plt.plot([coordsorig[l][0],coordsorig[s][0]],[coordsorig[l][1],coordsorig[s][1]],color=(0,0.0,0),linewidth=0.12*7*sfac)

    for k,o in nearbylst:
        h = False
        for p in lstall:
            if p[0] == k or p[1] == k:
                h = True
        if h == False:
            print('ERROR1 ',k)


    removedpointsall.append(lstv)
    print('\nremoving measured points')


    #closes point --
    if lstli == False:
        lstalli = []
        for a, b, i in coordsorig:
            for j in cnextnum[i]:
                ds = ((coordsorig[i][0] - coordsorig[j][0]) ** 2 + (coordsorig[i][1] - coordsorig[j][1]) ** 2) ** 0.5
                lstalli.append(((coordsorig[i][2]), (coordsorig[j][2]), ds))
    graphi = Graph(lstalli)  # redo

    if fix == True:
        coordclosest = ['a','b',routesfixc[cx-1]]
    if fix == False:
        coordclosest, dsmin = closest2(removedpointsall[-1], graphi)
    coordclosestroute = graphi.dijkstra(coordsorig[0][2], coordclosest[2])
    closestroute.append(coordclosestroute)

    closestpoint.append(coordclosest[2])
    #--


    removelst = []
    for l in coords:
        #print('test ',l)
        if l[2] in lstv:
            removelst.append(l)
    for i in removelst:
        print('removed ',i[2])
        coords.remove(i)
    print('')
    print(len(coords),' coords left')
    print('\n------ re-run ------')


    #run = False

#PLOT

#plot max radius
circle = plt.Circle((0,0),radm+100,color='r',fill=False,lw=0.5)
plt.gcf().gca().add_artist(circle)


#plot all  nodes
plotnodes = False
if plotnodes == True:
    plt.scatter([item[0] for item in coords],[item[1] for item in coords],s=12)


#plot runways
plt.plot([item[0] for item in rwy36Lcoords],[item[1] for item in rwy36Lcoords],'r')
plt.plot([item[0] for item in rwy24coords],[item[1] for item in rwy24coords],'r')
plt.plot([rwy36Lcoords[-1][0],rwy36Lcoords[0][0]],[rwy36Lcoords[-1][1],rwy36Lcoords[0][1]],'r')
plt.plot([rwy24coords[-1][0],rwy24coords[0][0]],[rwy24coords[-1][1],rwy24coords[0][1]],'r')


if f550 == True:
    removedpointsall = removedpointsall[:-1]

print('\n\nTotal flights: ',len(removedpointsall),'\n\n')
l1 = len(removedpointsall)
for c,rp in enumerate(removedpointsall):
    #cc = (c/l1,1-c/l1,abs(0.5-c/l1))
    cc = (random.randint(0,256)/256,random.randint(0,256)/256,random.randint(0,256)/256)#random colour per segment
    plt.scatter([coordsorig[item][0] for item in rp],[coordsorig[item][1] for item in rp],color=cc,s=60*sfac)#plot nodes per set
    plt.plot([coordsorig[item][0] for item in furthestroute[c]],[coordsorig[item][1] for item in furthestroute[c]],color=cc,linewidth=10*sfac,linestyle=':')#plot longest route
    plt.scatter([coordsorig[furthestpoint[c]][0]],[coordsorig[furthestpoint[c]][1]], color=cc, s=140*sfac)#plot longest route end node

    plt.plot([coordsorig[item][0] for item in closestroute[c]],[coordsorig[item][1] for item in closestroute[c]],color=cc,linewidth=7*sfac)#plot shortest route
    plt.scatter([coordsorig[closestpoint[c]][0]], [coordsorig[closestpoint[c]][1]], color=(0,0,0),s=190*sfac)  # plot shortest route end node
    plt.scatter([coordsorig[closestpoint[c]][0]], [coordsorig[closestpoint[c]][1]], color=cc, s=130 * sfac)  # plot shortest route end node

#plot ground station
plt.scatter([coordsorig[0][0]],[coordsorig[0][1]],color=(204/256,0,0),s=30)
plt.scatter([coordsorig[0][0]],[coordsorig[0][1]],color=(255/256,128/256,0),s=18)


#POINT CHECKER
checkval = int(12)
# plt.scatter([item[0] for item in cnext[checkval]],[item[1] for item in cnext[checkval]])
# plt.scatter([coordsorig[checkval][0]],[coordsorig[checkval][1]])



plt.axis('equal')
axes = plt.gca()
axes.set_xlim([-6000,6000])
axes.set_ylim([-6000,6000])
plt.show()
print("--- %s seconds ---" % (time.time() - start_time))

#SD
#plt.savefig('2_PP_SegmentedXXXX.png',dpi=600)

#SD2
#plt.savefig('2_PP_Segmented_V2_550_2.png',dpi=600)

#HD
#plt.savefig('2_PP_SegmentedXXXX.png',dpi=1000)
