import numpy as np
import cv2, PIL
from cv2 import aruco
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from PIL import Image

image_frame = cv2.imread("ArucoCode.JPG")
plt.figure()
plt.imshow(image_frame)
plt.axis('off')
#plt.show()
img = Image.open("ArucoCode.JPG")
width, height    = img.size
centre_x = width/2
centre_y = height/2

centre = [centre_x, centre_y]

gray = cv2.cvtColor(image_frame, cv2.COLOR_BGR2GRAY)
aruco_dict = aruco.Dictionary_get(aruco.DICT_6X6_250)
parameters =  aruco.DetectorParameters_create()
corners, marker, rejectedImgPoints = aruco.detectMarkers(gray, aruco_dict, parameters=parameters)
frame_markers = aruco.drawDetectedMarkers(image_frame.copy(), corners, marker)

plt.figure()
plt.imshow(frame_markers)
for i in range(len(marker)):
    c = corners[i][0]
    plt.plot([c[:, 0].mean()], [c[:, 1].mean()], "o", label = "id={0}".format('marker centre'))
plt.legend()

#plt.show()
def quad_area(data):
    l = data.shape[0]//2
    corners = data[["c1", "c2", "c3", "c4"]].values.reshape(l, 2,4)
    c1 = corners[:, :, 0]
    c2 = corners[:, :, 1]
    c3 = corners[:, :, 2]
    c4 = corners[:, :, 3]
    e1 = c2-c1
    e2 = c3-c2
    e3 = c4-c3
    e4 = c1-c4
    a = -.5 * (np.cross(-e1, e2, axis = 1) + np.cross(-e3, e4, axis = 1))
    return a, (e1[0][0] -e3[0][0] + e2[0][1] - e4[0][1])/4

corners2 = np.array([c[0] for c in corners])

data = pd.DataFrame({"x": corners2[:,:,0].flatten(), "y": corners2[:,:,1].flatten()},
                   index = pd.MultiIndex.from_product(
                           [marker.flatten(), ["c{0}".format(i )for i in np.arange(4)+1]],
                       names = ["marker", ""] ))



data = data.unstack().swaplevel(0, 1, axis = 1).stack()
data["m1"] = data[["c1", "c2"]].mean(axis = 1)
data["m2"] = data[["c2", "c3"]].mean(axis = 1)
data["m3"] = data[["c3", "c4"]].mean(axis = 1)
data["m4"] = data[["c4", "c1"]].mean(axis = 1)
data["o"] = data[["m1", "m2", "m3", "m4"]].mean(axis = 1)
data


corners2 = np.array([c[0] for c in corners])

data = pd.DataFrame({"x": corners2[:,:,0].flatten(), "y": corners2[:,:,1].flatten()},
                   index = pd.MultiIndex.from_product(
                           [marker.flatten(), ["c{0}".format(i )for i in np.arange(4)+1]],
                       names = ["marker", ""] ))

data = data.unstack().swaplevel(0, 1, axis = 1).stack()
data["m1"] = data[["c1", "c2"]].mean(axis = 1)
data["m2"] = data[["c2", "c3"]].mean(axis = 1)
data["m3"] = data[["c3", "c4"]].mean(axis = 1)
data["m4"] = data[["c4", "c1"]].mean(axis = 1)
data["o"] = data[["m1", "m2", "m3", "m4"]].mean(axis = 1)
#print('data:{}'.format(data))


Area, marker_length = quad_area(data)
marker_size = 9.74/100
distance_per_pixel = marker_size/marker_length
distance_to_centre_px = [(data['o'][0] -  centre[0]), (centre[1] - data['o'][1][1])]
total_distance_y  = float(distance_per_pixel)*float(distance_to_centre_px[0])
total_distance_x  = float(distance_per_pixel)*float(distance_to_centre_px[1])
print('The UAS has to move {} [m] in the x-direction and {} [m] in the y-direction'.format(total_distance_x, total_distance_y))
print(int(centre[1]))

#img = cv2.imread("ArucoCode.JPG")
cv2.arrowedLine(frame_markers, (int(centre[0]), int(centre[1])), (int(data['o'][0]),int(centre[1])), (255,0,0), 2)
cv2.arrowedLine(frame_markers, (int(data['o'][0]),int(centre[1])), (int(data['o'][0]), int(data['o'][1][1])), (255,0,0), 2)
plt.imshow(frame_markers)
plt.axis('off')
plt.show()