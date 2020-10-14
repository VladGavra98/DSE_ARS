import numpy as np
import cv2, PIL
from cv2 import aruco
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
#%matplotlib nbagg

aruco_dict = aruco.Dictionary_get(aruco.DICT_6X6_250)

fig = plt.figure()

img = aruco.drawMarker(aruco_dict,1, 700)
plt.imshow(img, cmap = mpl.cm.gray, interpolation = "nearest")
plt.axis("off")
plt.show()

