import numpy as np 
import math as m 
import matplotlib.pyplot as plt

def linear_accelerations(phi, theta, psi, T, g, m):
    sph = np.sin(phi)
    st  = np.sin(theta)
    sps = np.sin(psi)

    cph = np.cos(phi)
    ct = np.cos(theta)#
    cps = np.cos(psi)
    accelerations = np.add(g*np.array([0,0,1]), T/m*np.array([(cps*st*cph+sps*sph), (cps*sph-sps*st*cph), (-ct*cph)]))
    
    return accelerations

T = 1
m = 1
g = 9.81

phi = np.linspace(-np.pi/4,np.pi/4,100)
theta = np.linspace(-np.pi/4,np.pi/4,100)
psi = np.linspace(-np.pi/4,np.pi/4,100)

x_acc1 = np.array([])
y_acc1 = np.array([])
z_acc1 = np.array([])

for p in phi:
    t = 0
    ps = 0
    x_acc1 =np.append(x_acc1, linear_accelerations(p, t, ps, T, g, m)[0])
    y_acc1 = np.append(y_acc1, linear_accelerations(p, t, ps, T, g, m)[1])
    z_acc1 = np.append(z_acc1, linear_accelerations(p, t, ps, T, g, m)[2])

x_acc2 = np.array([])
y_acc2 = np.array([])
z_acc2 = np.array([])

for th in theta:
    p = 0
    ps = 0
    x_acc2 =np.append(x_acc2, linear_accelerations(p, th, ps, T, g, m)[0])
    y_acc2 = np.append(y_acc2, linear_accelerations(p, th, ps, T, g, m)[1])
    z_acc2 = np.append(z_acc2, linear_accelerations(p, th, ps, T, g, m)[2])

x_acc3 = np.array([])
y_acc3 = np.array([])
z_acc3 = np.array([])

for ps in psi:
    p = 0
    th = 0
    x_acc3 =np.append(x_acc3, linear_accelerations(p, th, ps, T, g, m)[0])
    y_acc3 = np.append(y_acc3, linear_accelerations(p, th, ps, T, g, m)[1])
    z_acc3 = np.append(z_acc3, linear_accelerations(p, th, ps, T, g, m)[2])


# plt.plot(phi, x_acc1, label = 'x_acc')
# plt.plot(phi, y_acc1, label = 'y_acc')
# plt.plot(phi, z_acc1, label = 'z_acc')
# plt.title('Relationship between linear accelerations and roll angle')
# plt.xlabel('roll angle')
# plt.ylabel('linear acceleration')
# plt.legend()
# plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(131)
ax1.plot(phi, x_acc1, label = 'x_acc')
ax1.plot(phi, y_acc1, label = 'y_acc')
ax1.plot(phi, z_acc1, label = 'z_acc')
ax1.set_xlabel('roll angle')
ax1.set_ylabel('linear acceleration')
ax1.legend()
ax2 = fig.add_subplot(132)
ax2.plot(theta, x_acc2, label = 'x_acc')
ax2.plot(theta, y_acc2, label = 'y_acc')
ax2.plot(theta, z_acc2, label = 'z_acc')
ax2.set_xlabel('pitch angle')
ax2.set_ylabel('linear acceleration')
ax2.legend()
ax3 = fig.add_subplot(133)
ax3.plot(psi, x_acc3, label = 'x_acc')
ax3.plot(psi, y_acc3, label = 'y_acc')
ax3.plot(psi, z_acc3, label = 'z_acc')
ax3.set_xlabel('yaw angle')
ax3.set_ylabel('linear acceleration')
ax3.legend()
fig.suptitle('Relationship between linear accelerations and Euler angles')
plt.show()

