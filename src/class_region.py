import numpy as np
from scipy.interpolate import NearestNDInterpolator
from class_sphere import Sphere
from class_cone import Cone

class Region:
    def __init__(self, Ntheta, Nphi, Nr, h_min, h_max, dr, r_pos, open_angle=30., region_in_ydirection=False, region_in_xdirection=True, upperhalf=True, cone_angle=30.):
        print('------------------Processing Region-------------------')
        assert not (region_in_ydirection==True and region_in_xdirection==True), 'Region object cannot be orientated in both x and y direction!'
        if open_angle < 90.:
            self.cone = Cone(Ntheta, Nr, h_min, h_max, dr, r_pos, open_angle, region_in_ydirection, region_in_xdirection, upperhalf)
            self.top = Sphere(Ntheta, Nphi, h_min, h_max, dr, True, r_pos, open_angle, region_in_ydirection, region_in_xdirection, upperhalf)
            self.bottom = Sphere(Ntheta, Nphi, h_min, h_max, dr, False, r_pos, open_angle, region_in_ydirection, region_in_xdirection, upperhalf)
        else:
            print('Cone object needs to have an open angle smaller than 90 deg. Initiate cone object with open angle = {} deg instead.'.format(cone_angle))
            self.cone = Cone(Ntheta, Nr, h_min, h_max, dr, r_pos, cone_angle, region_in_ydirection, region_in_xdirection, upperhalf)
            self.top = Sphere(Ntheta, Nphi, h_min, h_max, dr, True, r_pos, open_angle, region_in_ydirection, region_in_xdirection, upperhalf, cone_angle=cone_angle)
            self.bottom = Sphere(Ntheta, Nphi, h_min, h_max, dr, False, r_pos, open_angle, region_in_ydirection, region_in_xdirection, upperhalf, cone_angle=cone_angle)
        self.region_in_ydirection = region_in_ydirection
        self.region_in_xdirection = region_in_xdirection

    def get_region_dV_volume_exact(self):
        return self.top.dV_volume_exact, self.bottom.dV_volume_exact, self.cone.dV_volume_exact

    def get_region_dV_volume_approx(self, val, vol, count_jetr=False, verbose=False):
        self.top.get_surface_data(val, volume=vol, count_jetr=count_jetr)
        self.bottom.get_surface_data(val, volume=vol, count_jetr=count_jetr)
        self.cone.get_surface_data(val, volume=vol, count_jetr=count_jetr)

        top_vol = np.sum(self.top.volume)
        bot_vol = np.sum(self.bottom.volume)
        con_vol = np.sum(self.cone.volume)

        top_val = np.sum(self.top.volume * self.top.val)
        bot_val = np.sum(self.bottom.volume * self.bottom.val)
        con_val = np.sum(self.cone.volume * self.cone.val)

        if count_jetr:
            print('top count jetr:', self.top.counted_jetr)
            print('bot count jetr:', self.bottom.counted_jetr)
            print('con count jetr:', self.cone.counted_jetr)
            print('tot count jetr:', self.top.counted_jetr-self.bottom.counted_jetr+self.cone.counted_jetr)

        if verbose:
            print('top vol exact:', self.top.dV_volume_exact)
            print('top vol approx:', top_vol)
            print('bot vol exact:', self.bottom.dV_volume_exact)
            print('bot vol approx:', bot_vol)
            print('con vol exact:', self.cone.dV_volume_exact)
            print('con vol approx:', con_vol)

        # recalculate self.val to its surface value
        self.top.get_surface_data(val, volume=None)
        self.bottom.get_surface_data(val, volume=None)
        self.cone.get_surface_data(val, volume=None)

        return top_vol, bot_vol, con_vol, top_val, bot_val, con_val

    def get_region_surface_flux(self, value):
        self.top.surface_integral(value)
        self.flux_top = self.top.flux

        self.bottom.surface_integral(value)
        self.flux_bottom = self.bottom.flux

        self.cone.surface_integral(value)
        self.flux_cone = self.cone.flux
        
        self.total_flux = self.flux_top + self.flux_bottom + self.flux_cone

    def calc_Q_adv_conv_region(self, n, kT, vel):
        # mean value over surface:
        top_avg_n = self.top.surface_integral(n, avg=True)
        top_avg_kT = self.top.surface_integral(kT, avg=True)
        top_avg_vel_perp = self.top.surface_integral(vel, perp=True, avg=True)
        top_n = self.top.surface_integral(n, return_field=True)
        top_kT = self.top.surface_integral(kT, return_field=True)
        top_vel_perp = self.top.surface_integral(vel, perp=True, avg=False, return_field=True)

        bot_avg_n = self.bottom.surface_integral(n, avg=True)
        bot_avg_kT = self.bottom.surface_integral(kT, avg=True)
        bot_avg_vel_perp = self.bottom.surface_integral(vel, perp=True, avg=True)
        bot_n = self.bottom.surface_integral(n, return_field=True)
        bot_kT = self.bottom.surface_integral(kT, return_field=True)
        bot_vel_perp = self.bottom.surface_integral(vel, perp=True, avg=False, return_field=True)

        con_avg_n = self.cone.surface_integral(n, avg=True)
        con_avg_kT = self.cone.surface_integral(kT, avg=True)
        con_avg_vel_perp = self.cone.surface_integral(vel, perp=True, avg=True)
        con_n = self.cone.surface_integral(n, return_field=True)
        con_kT = self.cone.surface_integral(kT, return_field=True)
        con_vel_perp = self.cone.surface_integral(vel, perp=True, avg=False, return_field=True)

        # \delta value = value - mean(value)
        top_del_n = top_n - top_avg_n
        top_del_kT = top_kT - top_avg_kT
        top_del_vel_perp = top_vel_perp - top_avg_vel_perp
        #print(top_del_n.shape, top_del_vel_perp.shape, (top_del_n*top_del_vel_perp).shape)

        bot_del_n = bot_n - bot_avg_n
        bot_del_kT = bot_kT - bot_avg_kT
        bot_del_vel_perp = bot_vel_perp - bot_avg_vel_perp

        con_del_n = con_n - con_avg_n
        con_del_kT = con_kT - con_avg_kT
        con_del_vel_perp = con_vel_perp - con_avg_vel_perp

        # equation (11) from Yang+Reynolds2016:
        top_Q_adv  = 1.5 * (top_avg_n*top_avg_kT*top_avg_vel_perp + top_avg_kT*self.top.surface_integral(top_del_n*top_del_vel_perp, perp=False, avg=True, ignore_interpol=True))
        bot_Q_adv  = 1.5 * (bot_avg_n*bot_avg_kT*bot_avg_vel_perp + bot_avg_kT*self.bottom.surface_integral(bot_del_n*bot_del_vel_perp, perp=False, avg=True, ignore_interpol=True))
        con_Q_adv  = 1.5 * (con_avg_n*con_avg_kT*con_avg_vel_perp + con_avg_kT*self.cone.surface_integral(con_del_n*con_del_vel_perp, perp=False, avg=True, ignore_interpol=True))

        top_Q_conv = 1.5 * (top_avg_n*self.top.surface_integral(top_del_kT*top_del_vel_perp, perp=False, avg=True, ignore_interpol=True) + top_avg_vel_perp*self.top.surface_integral(top_del_n*top_del_kT, perp=False, avg=True, ignore_interpol=True) + self.top.surface_integral(top_del_n*top_del_kT*top_del_vel_perp, perp=False, avg=True, ignore_interpol=True))
        bot_Q_conv = 1.5 * (bot_avg_n*self.bottom.surface_integral(bot_del_kT*bot_del_vel_perp, perp=False, avg=True, ignore_interpol=True) + bot_avg_vel_perp*self.bottom.surface_integral(bot_del_n*bot_del_kT, perp=False, avg=True, ignore_interpol=True) + self.bottom.surface_integral(bot_del_n*bot_del_kT*bot_del_vel_perp, perp=False, avg=True, ignore_interpol=True))
        con_Q_conv = 1.5 * (con_avg_n*self.cone.surface_integral(con_del_kT*con_del_vel_perp, perp=False, avg=True, ignore_interpol=True) + con_avg_vel_perp*self.cone.surface_integral(con_del_n*con_del_kT, perp=False, avg=True, ignore_interpol=True) + self.cone.surface_integral(con_del_n*con_del_kT*con_del_vel_perp, perp=False, avg=True, ignore_interpol=True))

        # equation (10) from Yang+Reynolds2016:
        top_flux_L_adv = -top_Q_adv*self.top.area
        bot_flux_L_adv = -bot_Q_adv*self.bottom.area
        con_flux_L_adv = -con_Q_adv*self.cone.area

        top_flux_L_conv = -top_Q_conv*self.top.area
        bot_flux_L_conv = -bot_Q_conv*self.bottom.area
        con_flux_L_conv = -con_Q_conv*self.cone.area

        region_flux_L_adv = top_flux_L_adv + bot_flux_L_adv + con_flux_L_adv
        region_flux_L_conv = top_flux_L_conv + bot_flux_L_conv + con_flux_L_conv

        #return region_flux_L_adv, region_flux_L_conv
        return top_flux_L_adv, bot_flux_L_adv, con_flux_L_adv, top_flux_L_conv, bot_flux_L_conv, con_flux_L_conv

    def plot_xxyyzz_grid(self, snapn, bin, scale_fac):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(num=1)
        plt.clf()
        ax1 = fig.add_subplot(111, projection='3d')
        #ax1.view_init(elev=24., azim=-135.)
        npz_top = np.load('./npz/sphere_top_gridxxyyzz.npz')
        npz_bot = np.load('./npz/sphere_bot_gridxxyyzz.npz')
        npz_con = np.load('./npz/cone_gridxxyyzz.npz')

        xxt = npz_top['data_grid_xx']/scale_fac
        yyt = npz_top['data_grid_yy']/scale_fac
        zzt = npz_top['data_grid_zz']/scale_fac

        xxb = npz_bot['data_grid_xx']/scale_fac
        yyb = npz_bot['data_grid_yy']/scale_fac
        zzb = npz_bot['data_grid_zz']/scale_fac

        xxc = npz_con['data_grid_xx']/scale_fac
        yyc = npz_con['data_grid_yy']/scale_fac
        zzc = npz_con['data_grid_zz']/scale_fac

        if self.region_in_ydirection:
            ax1.scatter(xxc, zzc, yyc, s=0.2, alpha=0.5, color='r')
            ax1.scatter(xxb, zzb, yyb, s=0.2, alpha=0.5, color='g')
            ax1.scatter(xxt, zzt, yyt, s=0.2, alpha=0.5, color='b')
            ax1.set_xlabel(r'$xx$')
            ax1.set_ylabel(r'$zz$')
            ax1.set_zlabel(r'$yy$')
        elif self.region_in_xdirection:
            ax1.scatter(yyc, zzc, xxc, s=0.2, alpha=0.5, color='r')
            ax1.scatter(yyb, zzb, xxb, s=0.2, alpha=0.5, color='g')
            ax1.scatter(yyt, zzt, xxt, s=0.2, alpha=0.5, color='b')
            ax1.set_xlabel(r'$yy$')
            ax1.set_ylabel(r'$zz$')
            ax1.set_zlabel(r'$xx$')
        else:
            ax1.scatter(xxc, yyc, zzc, s=0.2, alpha=0.5, color='r')
            ax1.scatter(xxb, yyb, zzb, s=0.2, alpha=0.5, color='g')
            ax1.scatter(xxt, yyt, zzt, s=0.2, alpha=0.5, color='b')
            ax1.set_xlabel(r'$xx$')
            ax1.set_ylabel(r'$yy$')
            ax1.set_zlabel(r'$zz$')

        fig.tight_layout()
        plt.savefig('./images/3d_scatter_region_all_xxyyzz_n{:02d}_b{}.png'.format(snapn, bin), dpi=300)
        plt.close('all')
        print('Plotted xxyyzz grid.')

    def plot_xyz_grid(self, snapn, bin, scale_fac):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(num=1)
        plt.clf()
        ax1 = fig.add_subplot(111, projection='3d')
        #ax1.view_init(elev=24., azim=-135.)
        npz_top = np.load('./npz/sphere_top_gridxyz.npz')
        npz_bot = np.load('./npz/sphere_bot_gridxyz.npz')
        npz_con = np.load('./npz/cone_gridxyz.npz')

        xt = npz_top['data_grid_x']/scale_fac
        yt = npz_top['data_grid_y']/scale_fac
        zt = npz_top['data_grid_z']/scale_fac

        xb = npz_bot['data_grid_x']/scale_fac
        yb = npz_bot['data_grid_y']/scale_fac
        zb = npz_bot['data_grid_z']/scale_fac

        xc = npz_con['data_grid_x']/scale_fac
        yc = npz_con['data_grid_y']/scale_fac
        zc = npz_con['data_grid_z']/scale_fac

        if self.region_in_ydirection:
            ax1.scatter(xc, zc, yc, s=0.2, alpha=0.5, color='r')
            ax1.scatter(xb, zb, yb, s=0.2, alpha=0.5, color='g')
            ax1.scatter(xt, zt, yt, s=0.2, alpha=0.5, color='b')
            ax1.set_xlabel(r'$x$')
            ax1.set_ylabel(r'$z$')
            ax1.set_zlabel(r'$y$')
        elif self.region_in_xdirection:
            ax1.scatter(yc, zc, xc, s=0.2, alpha=0.5, color='r')
            ax1.scatter(yb, zb, xb, s=0.2, alpha=0.5, color='g')
            ax1.scatter(yt, zt, xt, s=0.2, alpha=0.5, color='b')
            ax1.set_xlabel(r'$y$')
            ax1.set_ylabel(r'$z$')
            ax1.set_zlabel(r'$x$')
        else:
            ax1.scatter(xc, yc, zc, s=0.2, alpha=0.5, color='r')
            ax1.scatter(xb, yb, zb, s=0.2, alpha=0.5, color='g')
            ax1.scatter(xt, yt, zt, s=0.2, alpha=0.5, color='b')
            ax1.set_xlabel(r'$x$')
            ax1.set_ylabel(r'$y$')
            ax1.set_zlabel(r'$z$')

        fig.tight_layout()
        plt.savefig('./images/3d_scatter_region_all_xyz_n{:02d}_b{}.png'.format(snapn, bin), dpi=300)
        plt.close('all')
        print('Plotted xyz grid.')