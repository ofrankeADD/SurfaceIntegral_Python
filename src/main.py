'''
author: Oliver Franke
email: ofranke@posteo.de
github: ofrankeADD
latest update: 25.11.2021
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from mpl_toolkits.mplot3d import Axes3D

def load_snap(filename):
    ptName = 'PartType0'
    fields = {}
    with h5py.File(filename, 'r') as f:
        for k in f[ptName].keys():
            fields[k] = f[ptName][k][()]
            #print(k)
        #time = f['Header'].attrs['Time']
        #for a in f['Header'].attrs.keys():
        #    print(a)
        #h67 = f['Header'].attrs['HubbeParam']
    return fields

def get_npz_data(use_h5py, rmax, snapnum, snappath):
    try:
        if use_h5py:
            npz_data = np.load('./dS_h5py.npz')
        else:
            npz_data = np.load('./dS_snap.npz')
        r = npz_data['data_r']
        r_pos = npz_data['data_rpos']
        rho = npz_data['data_rho']
        vel = npz_data['data_vel']
    except:
        if use_h5py:
            import h5py
            fields = load_snap(snappath+'/snap_{:03d}.hdf5'.format(snapnum))
            rho = np.array(fields['Density'][:])
            vel = np.array(fields['Velocity'][:])
            r_pos = fields['Coordinates'][:]/0.67 - np.array([0.5/0.67]*3)
            r = np.sqrt(np.sum(r_pos**2, axis=1))
            r = np.array(r*1e3)
            r_pos = np.array(r_pos*1e3)
            r_pos = r_pos[(r<rmax)]
            rho = rho[(r<rmax)]
            vel = vel[(r<rmax)]
            r = r[(r<rmax)]
            np.savez('./dS_h5py.npz', data_r=r, data_rpos=r_pos, data_rho=rho, data_vel=vel)
        else:
            import gadget
            snap = gadget.gadget_readsnap(snapnum, snapbase='snap_', snappath=snappath)
            rho = snap.data['rho']
            vel = snap.data['vel'][snap.data['type']==0]
            r_pos = snap.data['pos']-snap.center
            r = np.sqrt(np.sum(r_pos**2, axis=1))
            r = r[snap.data['type']==0]
            r *= 1e3
            r_pos = r_pos[snap.data['type']==0]
            r_pos *= 1e3
            r_pos = r_pos[(r<rmax)]
            rho = rho[(r<rmax)]
            vel = vel[(r<rmax)]
            r = r[(r<rmax)]
            np.savez('./dS_snap.npz', data_r=r, data_rpos=r_pos, data_rho=rho, data_vel=vel)
    return r, r_pos, rho, vel

'''

'''
if __name__ == "__main__":

# if False use gadget:
    use_h5py = False
# limiting the size of the npz file:
    rmax = 100 # kpc
# which snapshot to load:
    snapnum = 0
    snappath = './output_brag_beta20/'
    r, r_pos, rho, vel = get_npz_data(use_h5py, rmax, snapnum, snappath)

# defining a thin shell boundary at r=50 kpc with width dr=5 kpc:
    dr = 5. # kpc
    rbin = 50. # kpc
    shell_rho   =   rho[(rbin-dr/2. < r) & (r < rbin+dr/2.)]
    shell_vel   =   vel[(rbin-dr/2. < r) & (r < rbin+dr/2.)]
    shell_r     =     r[(rbin-dr/2. < r) & (r < rbin+dr/2.)]
    shell_pos   = r_pos[(rbin-dr/2. < r) & (r < rbin+dr/2.)]
    shell_vr = (np.sum(shell_vel*shell_pos, axis=1))/shell_r
    shell_v = np.sqrt(np.sum(shell_vel**2, axis=1))
    print(shell_r.shape, shell_pos.shape, shell_rho.shape, shell_vr.shape)
    print('xyz voronoi rho: mean:', np.mean(shell_rho))
    print('shell_vr:', np.min(shell_vr), np.max(shell_vr))

# getting the xyz coordinates and transforming them into spherical coordinates theta and phi:
# axes orientation according to https://en.wikipedia.org/wiki/Spherical_coordinate_system
    shell_x = shell_pos[:,2]
    shell_y = shell_pos[:,0]
    shell_z = shell_pos[:,1]
    shell_phi = np.arctan2(shell_y, shell_x)
    shell_theta = np.arccos(shell_z/shell_r)
    print(shell_x.shape, shell_theta.shape, shell_phi.shape)

# initiating a grid for the spherical interpolation, where theta in [0, pi] and phi in [0, 2*pi] according to wikipedia
# dtheta = pi/N and dphi = 2*pi/N
# radius is just a single value at rbin
    Nr = 1
    Nphi = 100
    Ntheta = 100
    #sph_r = np.array([rbin]*Nr)
    sph_r = rbin
    sph_phi = np.linspace(0., 2.*np.pi, Nphi, endpoint=True)
    sph_theta_incorrect_old = np.linspace(0., np.pi, Ntheta, endpoint=True)
    #sph_theta = np.linspace(0., np.pi, Ntheta, endpoint=True)
# why the above doesn't yield a uniform distribution of theta see:
# https://corysimon.github.io/articles/uniformdistn-on-sphere/
# https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
    correct_uniform_theta = np.linspace(0., 1., Ntheta, endpoint=True)
    sph_theta = np.arccos(1. - 2.*correct_uniform_theta)

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
    v_x = xx**2
    v_y = yy**2
    v_z = zz
    dS_vr = v_x*dS_x + v_y*dS_y + v_z*dS_z

    surf_int_approx = np.sum(np.sum(dS_vr))
    surf_int_exact = 8./3.*np.pi*sph_r**4
    #surf_int_exact = 4.*np.pi*sph_r**4

    print('Nr: {}, Ntheta: {}, Nphi: {}'.format(Nr, Ntheta, Nphi))
    print('Exact  surface integral: {:.1f}'.format(surf_int_exact))
    print('Approx surface integral: {:.1f}'.format(surf_int_approx))
    print('Difference: {:.1f}'.format(surf_int_exact - surf_int_approx))
    print('absolute error: {:.1f}%\n'.format(100.*np.abs(surf_int_exact - surf_int_approx)/np.maximum(surf_int_exact, surf_int_approx)))

# \int v\cdot r dS from actual velocity field:
    interp_vel = NearestNDInterpolator(list(zip(shell_x, shell_y, shell_z)), shell_vel)
    vel_interp = interp_vel(xx, yy, zz)
    print(vel_interp.shape)

    vr_interp1 = n_vec[:,:,0]*vel_interp[:,:,0] + n_vec[:,:,1]*vel_interp[:,:,1] + n_vec[:,:,2]*vel_interp[:,:,2]
    print(vr_interp1.shape)

    surf_int_vr = np.sum(np.sum(vr_interp1*dS))
    print('vr_interp1.shape:', vr_interp1.shape)
    print('surf_int_vr:', surf_int_vr)





#############################################################################################################

    plot_vr_interp = False
    if plot_vr_interp:
        #v_x = vr*np.sin(theta)*np.cos(phi)
        #v_y = vr*np.sin(theta)*np.sin(phi)
        #v_z = vr*np.cos(theta)

        fig = plt.figure(num=1)
        plt.clf()
        #fig, ax = plt.subplots(num=1, ncols=3)

        #vmin0 = np.min(v_xyz_r)
        #vmax0 = np.max(v_xyz_r*n_vec[:,:,0])
        #print(vmin0, vmax0)
        #vmin1 = np.min(np.sqrt(vr_xx**2+vr_yy**2+vr_zz**2))
        #vmax1 = np.max(np.sqrt(vr_xx**2+vr_yy**2+vr_zz**2))
        #print(vmin1, vmax1)
        vmin2 = np.min(vr_interp1)
        vmax2 = np.max(vr_interp1)
        print(vmin2, vmax2)

        cmap = plt.cm.YlOrRd

        ax1 = fig.add_subplot(121, projection='3d')
        #im2 = ax[1].pcolormesh(xx, yy, np.sqrt(vr_xx**2+vr_yy**2+vr_zz**2), vmin=-200., vmax=200.)
        #ax1.view_init(elev=24., azim=-135.)
        ax1.scatter(xx, yy, zz, s=0.4, alpha=0.5)
        ax1.set_xlabel(r'$x$')
        ax1.set_ylabel(r'$y$')
        ax1.set_zlabel(r'$z$')

        ax2 = fig.add_subplot(122)
        im2 = ax2.pcolormesh(xx, yy, vr_interp1, vmin=-200., vmax=200.)#, cmap=cmap)
        ax2.set_xlabel(r'$x$')
        ax2.set_ylabel(r'$y$')

        #ax3 = fig.add_subplot(133, projection='3d')
        #im3 = ax3.pcolormesh(xx, yy, v_xyz_r*n_vec[:,:,0], vmin=-200., vmax=200.)
        ##im1 = ax[0].scatter(shell_r, shell_vr, s=0.8)
        #im3 = ax3.plot_surface(xx, yy, zz, rstride=1, cstride=1, facecolors=cmap(vr_interp1), shade=False)
        #im3 = ax3.plot_surface(xx, yy, zz, vmin=-200, vmax=200, rstride=1, cstride=1, color=cmap(vr_interp1))
        
        fig.tight_layout()
        plt.savefig('./pcolor_vr.png', dpi=300)
        plt.close('all')




    plot_G_profile = False
    if plot_G_profile:

        G_analyt = np.cos(tt) * np.sin(4.*pp)
        G_shell = np.cos(shell_theta) * np.sin(4.*shell_phi)
        print(G_shell.shape, G_shell)

        interp_G = NearestNDInterpolator(list(zip(shell_x, shell_y, shell_z)), G_shell)
        G_interp = interp_G(xx, yy, zz)

        #interp_G_ana = NearestNDInterpolator(list(zip(sph_theta, sph_phi)), G_analyt)
        #G_interp_ana = interp_G_ana(xx, yy)
        #interp_G = NearestNDInterpolator(list(zip(shell_theta, shell_phi)), G_shell)
        #G_interp = interp_G(tt, pp)
        print(G_interp.shape)

    # plot using pcolormesh:
        fig = plt.figure(num=1)
        plt.clf()
        fig, ax = plt.subplots(num=1, ncols=2)

        extent_tp = [sph_theta[0], sph_theta[-1], sph_phi[0], sph_phi[-1]]
        extent_xy = [-rbin, rbin, -rbin, rbin]
        print(extent_tp)
        vmin1 = np.min(G_interp)
        vmax1 = np.max(G_interp)
        #vmin2 = np.min(G_interp_ana)
        #vmax2 = np.max(G_interp_ana)
        print(vmin1, vmax1)
        #print(vmin2, vmax2)
        ax[0].autoscale(False)
        ax[0].set_xlim(extent_xy[0], extent_xy[1])
        ax[0].set_ylim(extent_xy[2], extent_xy[3])
        im0 = ax[0].pcolormesh(xx, yy, G_interp, vmin=-1., vmax=1.)
        cb0 = fig.colorbar(im0, ax=ax[0], orientation='horizontal', shrink=0.75)
        cb0.set_label('G_interp')
        #ax[0].set_xlabel(r'$\mathrm{theta}\ \theta$')
        #ax[0].set_ylabel(r'$\mathrm{phi}\ \phi$')
        ax[0].set_xlabel(r'$x$')
        ax[0].set_ylabel(r'$y$')

        ax[1].autoscale(False)
        ax[1].set_xlim(extent_xy[0], extent_xy[1])
        ax[1].set_ylim(extent_xy[2], extent_xy[3])
        im1 = ax[1].pcolormesh(xx, yy, G_analyt, vmin=-1., vmax=1.)
        cb1 = fig.colorbar(im1, ax=ax[1], orientation='horizontal', shrink=0.75)
        cb1.set_label('G_analyt')
        #ax[1].set_xlabel(r'$\mathrm{theta}\ \theta$')
        #ax[1].set_ylabel(r'$\mathrm{phi}\ \phi$')
        ax[1].set_xlabel(r'$x$')
        ax[1].set_ylabel(r'$y$')

    # plot 1D analytical profile instead:
        #ax[1].plot(G_analyt, color='k', linestyle='dotted')
        #ax[1].set_xlabel('')
        #ax[1].set_ylabel('G_analyt')

        fig.tight_layout()
        plt.savefig('./pcolor_G_xy.png', dpi=300)
        plt.close('all')


