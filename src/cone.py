import numpy as np
import matplotlib.pyplot as plt

"""
This python script sets up a grid on a truncated cone.

Relevant links are:

https://math.stackexchange.com/questions/855831/surface-integral-of-a-right-circular-cone

A useful test case would be reproducing the problem described here

https://math.stackexchange.com/questions/2257382/finding-the-unit-normal-to-a-cone

both numerically and analytically.

written by Thomas Berlok on 29/1/2021
"""

# Opening angle of cone
alpha = np.pi/6

# Minimum radius
r_min = 1

# Maximum radius
r_max = 2

# Corresponding heights above z plane
h_min = r_min/np.tan(alpha)
h_max = r_max/np.tan(alpha)

# Resolution for the grid
Nr = 400
Ntheta = 400
dr = (r_max - r_min)/Nr
dtheta = 2*np.pi/Ntheta

# Grid 1D arrays
r_vec = r_min + (np.arange(Nr) + 1/2) * dr
theta_vec  = (np.arange(Ntheta) + 1/2) * dtheta

# Grid 2d arrays
r_mat, theta_mat = np.meshgrid(r_vec, theta_vec)

# height 2d array
hh = r_mat/np.tan(alpha)

# Cartesian grid arrays
xx = r_mat*np.cos(theta_mat)
yy = r_mat*np.sin(theta_mat)
zz = hh

# Components of normal vector 2d array
n_x = -np.cos(theta_mat)/np.tan(alpha)
n_y = -np.sin(theta_mat)/np.tan(alpha)
n_z = np.ones_like(r_mat)

# Area element 2d arrays. \vec{dS} =  \vec{n} dx dy = \vec{n} r dr dθ
dS_x = n_x * r_mat * dr * dtheta
dS_y = n_y * r_mat * dr * dtheta
dS_z = n_z * r_mat * dr * dtheta

# Area of cone stub
dS = np.sqrt(dS_x**2 + dS_y**2 + dS_z**2)
area_numerical = np.sum(np.sum(dS))

# Exact result from here
# https://rechneronline.de/pi/truncated-cone.php
R = r_max
r = r_min
h = h_max - h_min
s = np.sqrt( (R-r)**2 + h**2 )
area_exact = (R+r)*np.pi*s
print('Nr: {}, Ntheta: {}'.format(Nr, Ntheta))
print('Exact cone area: {}'.format(area_exact))
print('Approximate cone area: {}'.format(area_numerical))
print('Difference: {}\n'.format(area_exact-area_numerical))

# A somewhat arbitrary velocity field for testing
v_x = xx**3
v_y = yy**3
v_z = zz

# v ⋅ dS
tmp = v_x*dS_x + v_y*dS_y + v_z*dS_z

surface_integral_numerical = np.sum(tmp[:])
surface_integral_exact = np.pi/(30*np.tan(alpha))*(20*(r_max**3 - r_min**3) - 9*(r_max**5 - r_min**5))

print('Nr: {}, Ntheta: {}'.format(Nr, Ntheta))
print('Exact surface integral: {}'.format(surface_integral_exact))
print('Approximate surface integral: {}'.format(surface_integral_numerical))
print('Difference: {}\n'.format(surface_integral_exact-surface_integral_numerical))
