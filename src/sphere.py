import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

# initiating a grid for the spherical interpolation, where theta in [0, pi] and phi in [0, 2*pi] according to wikipedia
# dtheta = pi/N and dphi = 2*pi/N
# radius is just a single value at rbin
    rbin = 50 # kpc
    Nr = 1
    Nphi = 200
    Ntheta = 200
    #sph_r = np.array([rbin]*Nr)
    sph_r = rbin
    sph_phi = np.linspace(0., 2.*np.pi, Nphi, endpoint=True)
    sph_theta_incorrect_old = np.linspace(0., np.pi, Ntheta, endpoint=True)
    sph_theta = np.linspace(0., np.pi, Ntheta, endpoint=True)
# why the above doesn't yield a uniform distribution of theta see:
# https://corysimon.github.io/articles/uniformdistn-on-sphere/
# https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
    correct_uniform_theta = np.linspace(0., 1., Ntheta, endpoint=True)
    #sph_theta = np.arccos(1. - 2.*correct_uniform_theta)

# compare analytical profile with interpolated one:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.NearestNDInterpolator.html
    tt, pp = np.meshgrid(sph_theta, sph_phi)
    xx = sph_r*np.sin(tt)*np.cos(pp)
    yy = sph_r*np.sin(tt)*np.sin(pp)
    zz = sph_r*np.cos(tt)

    n_vec = np.zeros((Ntheta,Nphi,3))
    n_vec[:,:,0] = xx
    n_vec[:,:,1] = yy
    n_vec[:,:,2] = zz
    n_vec /= sph_r

    # n_x = sph_r**2 * np.sin(sph_theta)**2 * np.cos(sph_phi)
    # n_y = sph_r**2 * np.sin(sph_theta)**2 * np.sin(sph_phi)
    # n_z = -sph_r**2 * np.cos(sph_theta) * np.sin(sph_theta)

# elementary surface area dS = r**2 * sin(theta) * dtheta * dphi
    dtheta  = np.pi/Ntheta
    dphi  = 2.*np.pi/Nphi
    dS = sph_r**2 * np.sin(tt) * dtheta * dphi
    dS_x = n_vec[:,:,0] * sph_r**2 * np.sin(tt) * dtheta * dphi
    dS_y = n_vec[:,:,1] * sph_r**2 * np.sin(tt) * dtheta * dphi
    dS_z = n_vec[:,:,2] * sph_r**2 * np.sin(tt) * dtheta * dphi

    # dS_x = n_x * sph_r**2 * np.sin(sph_theta) * dtheta * dphi
    # dS_y = n_y * sph_r**2 * np.sin(sph_theta) * dtheta * dphi
    # dS_z = n_z * sph_r**2 * np.sin(sph_theta) * dtheta * dphi

    dS_xyz = np.sqrt(dS_x**2 + dS_y**2 + dS_z**2)
    dS_xyz_area = np.sum(np.sum(dS_xyz))
    print('dS_xyz: mean:', np.mean(dS_xyz), dS_xyz.shape)
    print('dS_xyz: sum:', dS_xyz_area)
    print('dS: sum:', np.sum(dS))

# compare exact with approx surface area:
    dS_area_exact = 4.*np.pi*sph_r**2
    print('Nr: {}, Ntheta: {}, Nphi: {}'.format(Nr, Ntheta, Nphi))
    print('Exact  sphere area: {:.1f}'.format(dS_area_exact))
    print('Approx sphere area: {:.1f}'.format(dS_xyz_area))
    print('Difference: {:.1f}'.format(dS_area_exact - dS_xyz_area))
    print('absolute error: {:.1f}%\n'.format(100.*np.abs(dS_area_exact - dS_xyz_area)/np.maximum(dS_area_exact, dS_xyz_area)))

# arbitrary velocity field:
    v_x = xx**3
    v_y = yy**3
    v_z = zz
    dS_vr = v_x*dS_x + v_y*dS_y + v_z*dS_z
    # dS_vr = (v_x*n_vec[:,:,0] + v_y*n_vec[:,:,1] + v_z*n_vec[:,:,2]) * sph_r**2 * np.sin(tt) * dtheta * dphi

    surf_int_approx = np.sum(np.sum(dS_vr))
    surf_int_exact = 8./3.*np.pi*sph_r**4
    surf_int_exact = 8.*(np.pi*(4.*sph_r+3.)/(12.*sph_r))*sph_r**4
    surf_int_exact = 4./3.*np.pi*sph_r

    # The integral below is calculated using Mathematica as follows:
    # x = r Sin[\[Theta]] Cos[\[Phi]];
    # y = r Sin[\[Theta]] Sin[\[Phi]];
    # z = r Cos[\[Theta]];
    # n = {x, y, z}/r;
    # dS = r^2 Sin[\[Theta]] n ;
    # V = {x^3, y^3, z};
    # Integrate[Integrate[V.dS, {\[Phi], 0, 2 \[Pi]}], {\[Theta], 0, \[Pi]}]
    surf_int_exact = 4/15*np.pi*sph_r**3*(5 + 6*sph_r**2)

    # Note also that for
    # v_x = xx**2
    # v_y = yy**2
    # v_z = zz
    # i.e. with squares instead of cubes (that you use for some reason),
    # the exact result is
    # surf_int_exact = 4/3*np.pi*sph_r**3
    # This version converges but requires many more points!

    print('Nr: {}, Ntheta: {}, Nphi: {}'.format(Nr, Ntheta, Nphi))
    print('Exact  surface integral: {:.1f}'.format(surf_int_exact))
    print('Approx surface integral: {:.1f}'.format(surf_int_approx))
    print('Difference: {:.1f}'.format(surf_int_exact - surf_int_approx))
    print('absolute error: {:.1f}%\n'.format(100.*np.abs(surf_int_exact - surf_int_approx)/np.maximum(surf_int_exact, surf_int_approx)))



