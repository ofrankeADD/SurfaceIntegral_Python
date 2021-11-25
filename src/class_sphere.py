import numpy as np
from scipy.interpolate import NearestNDInterpolator
from class_surface import Surface

class Sphere(Surface):
    def __init__(self, Ntheta, Nphi, h_min, h_max, dr, n_top, r_pos, open_angle=30., sphere_in_ydirection=False, sphere_in_xdirection=True, upperhalf=True, cone_angle=30.):
        print('Initiating Sphere object.')
        self.type = 'Sphere'
        self.n_top = n_top
        self.open_angle = np.radians(open_angle)

        if not self.n_top:
            if self.open_angle < np.pi/2:
                self.r = h_min + h_min*(1. - np.cos(self.open_angle))
                self.h = h_min
            else:
                self.cone_angle = np.radians(cone_angle)
                self.r = h_min + h_min*(1. - np.cos(self.cone_angle))
                self.h = 0.
        else:
            if self.open_angle < np.pi/2:
                self.r = h_max + h_max*(1. - np.cos(self.open_angle))
                self.h = h_max
            else:
                self.cone_angle = np.radians(cone_angle)
                self.r = h_max + h_max*(1. - np.cos(self.cone_angle))
                self.h = 0.

        self.dr = dr
        self.flux = 0.
        self.Nr = 1
        self.Ntheta = Ntheta
        self.Nphi = Nphi
        self.r_pos = r_pos
        assert not (sphere_in_ydirection==True and sphere_in_xdirection==True), 'Sphere object cannot be orientated in both x and y direction!'
        self.sphere_in_ydirection = sphere_in_ydirection
        self.sphere_in_xdirection = sphere_in_xdirection
        self.upperhalf = upperhalf
        self.get_dS_area_exact()
        self.get_dV_volume_exaxt()

        self.create_grid()
        Surface.__init__(self)

    def get_dS_area_exact(self):
        if self.open_angle < np.pi/2.:
            self.dS_area_exact = 2.*np.pi*self.r**2 * (1. - np.cos(self.open_angle))
        else:
            self.dS_area_exact = 0.5 * 4.*np.pi * self.r**2

    def get_dV_volume_exaxt(self):
        if self.open_angle < np.pi/2.:
            self.dV_volume_exact = (np.pi/3.) * self.r**3 * (2. + np.cos(self.open_angle)) * (1. - np.cos(self.open_angle))**2
        else:
            self.dV_volume_exact = 0.5 * 4./3. * np.pi * self.r**3

    def iterate_surface_integral(self, value, surf_int_compare, stepsize=100, tol=0.01, perp=False):
        err = Surface.iterate_surface_integral(self, value, surf_int_compare, stepsize, tol)
        return err

    def iterate_grid_resolution(self, stepsize=100, tol=0.01):
        Surface.iterate_grid_resolution(self, stepsize, tol)

    def surface_integral(self, value, perp=False, avg=False, return_field=False, ignore_interpol=False, verbose=False):
        if verbose:
            print('Calculating Surface integral for {} object.'.format(self.type))

        if not ignore_interpol:
            self.get_surface_data(value)
            #self.flux = Surface.surface_integral(self)
            field = self.interpolate_to_grid()
        else:
            field = value

        is_vector_field = (len(field.shape) == 3 and field.shape[-1] == 3)

        if perp:
            assert is_vector_field, 'perp set to True but input is not a vector field'
            field = field[:,:,0]*self.n_vec[:,:,0] + field[:,:,1]*self.n_vec[:,:,1] + field[:,:,2]*self.n_vec[:,:,2]
            is_vector_field = len(field.shape) == 3 and field.shape[-1] == 3

        if return_field:
            return field

        dS_x = self.dS_vec[:,:,0]
        dS_y = self.dS_vec[:,:,1]
        dS_z = self.dS_vec[:,:,2]
        self.dS = np.sqrt(dS_x**2 + dS_y**2 + dS_z**2)
        self.area = np.sum(self.dS)

        if is_vector_field:
            vec_x = field[:,:,0]
            vec_y = field[:,:,1]
            vec_z = field[:,:,2]

            self.flux = np.sum(vec_x*dS_x + vec_y*dS_y + vec_z*dS_z)
        else:
            #self.flux = np.sum(dS_x*field + dS_y*field + dS_z*field)
            self.flux = np.sum(self.dS*field)

        if avg:
            if verbose:
                print('average surface flux: {}\n'.format(self.flux/self.area))
            return self.flux/self.area
        else:
            if verbose:
                print('integrated surface flux: {}\n'.format(self.flux))
            return self.flux

    def interpolate_to_grid(self, verbose=False):
        if verbose:
            print('Interpolating Surface data on grid for {} object.'.format(self.type))
        interp_v = NearestNDInterpolator(list(zip(self.x, self.y, self.z)), self.val)
        v_interp = interp_v(self.xx, self.yy, self.zz)

        return v_interp

    def create_grid(self, alternative_theta=False, verbose=False):
        if verbose:
            print('Creating grid for Sphere object.')
        # initiating a grid for the spherical interpolation, where theta in [0, open_angle] and phi in [0, 2*pi] according to wikipedia
        # radius is just a single value at r
        sph_phi = np.linspace(0., 2.*np.pi, self.Nphi, endpoint=False)
        if alternative_theta:
            sph_theta = np.linspace(-self.open_angle, self.open_angle, self.Ntheta, endpoint=True)
        else:
            sph_theta = np.linspace(0, self.open_angle, self.Ntheta, endpoint=False)

        tt, pp = np.meshgrid(sph_theta, sph_phi)
        self.xx = self.r*np.sin(tt)*np.cos(pp)
        self.yy = self.r*np.sin(tt)*np.sin(pp)
        self.zz = self.r*np.cos(tt)
        
        if not self.sphere_in_xdirection:
            if self.upperhalf:
                self.zz *= -1.

        n_vec = np.zeros((self.Nphi,self.Ntheta,3))
        n_vec[:,:,0] = self.xx
        n_vec[:,:,1] = self.yy
        n_vec[:,:,2] = self.zz
        n_vec /= self.r
        if not self.n_top:
            n_vec *= -1.
        self.n_vec = n_vec

        # elementary surface area d\vec{S} = \vec{n} * r**2 * sin(theta) * dtheta * dphi
        if alternative_theta:
            dtheta  = self.open_angle/(self.Ntheta-1.)
        else:
            dtheta  = self.open_angle/self.Ntheta

        dphi  = 2.*np.pi/self.Nphi
        dS_vec = np.zeros((self.Nphi,self.Ntheta,3))
        dS_vec[:,:,0] = self.n_vec[:,:,0] * self.r**2 * np.sin(tt) * dtheta * dphi
        dS_vec[:,:,1] = self.n_vec[:,:,1] * self.r**2 * np.sin(tt) * dtheta * dphi
        dS_vec[:,:,2] = self.n_vec[:,:,2] * self.r**2 * np.sin(tt) * dtheta * dphi
        self.dS_vec = dS_vec

        if self.sphere_in_ydirection:
            self.rotate_90deg_around_x(verbose=False)
        elif self.sphere_in_xdirection:
            self.rotate_90deg_around_y(verbose=False)

        if self.n_top:
            np.savez('./npz/sphere_top_gridxxyyzz.npz', data_grid_xx=self.xx, data_grid_yy=self.yy, data_grid_zz=self.zz)
        else:
            np.savez('./npz/sphere_bot_gridxxyyzz.npz', data_grid_xx=self.xx, data_grid_yy=self.yy, data_grid_zz=self.zz)

    def rotate_90deg_around_x(self, verbose=False):
        tmp_yy = np.copy(self.yy)
        self.yy = -np.copy(self.zz)
        self.zz = tmp_yy

        tmp_yy = np.copy(self.n_vec[:,:,1])
        self.n_vec[:,:,1] = -np.copy(self.n_vec[:,:,2])
        self.n_vec[:,:,2] = tmp_yy

        tmp_yy = np.copy(self.dS_vec[:,:,1])
        self.dS_vec[:,:,1] = -np.copy(self.dS_vec[:,:,2])
        self.dS_vec[:,:,2] = tmp_yy

        if verbose:
            print('Sphere object is now orientated in y direction.')

    def rotate_90deg_around_y(self, verbose=False):
        tmp_xx = np.copy(self.xx)
        self.xx = np.copy(self.zz)
        self.zz = -tmp_xx

        tmp_xx = np.copy(self.n_vec[:,:,0])
        self.n_vec[:,:,0] = np.copy(self.n_vec[:,:,2])
        self.n_vec[:,:,2] = -tmp_xx

        tmp_xx = np.copy(self.dS_vec[:,:,0])
        self.dS_vec[:,:,0] = np.copy(self.dS_vec[:,:,2])
        self.dS_vec[:,:,2] = -tmp_xx

        if verbose:
            print('Cone object is now orientated in x direction.')

    def get_surface_data(self, value, volume=None, count_jetr=False, verbose=False):
        if verbose:
            print('Initiating Surface data for {} object.'.format(self.type))

        r_len = np.sqrt(np.sum(self.r_pos**2, axis=1))

        # defining a thin shell boundary at r kpc with width dr kpc:
        sel_crit = ((self.r-self.dr/2. < r_len) & (r_len < self.r+self.dr/2.))
        shell_value = value[sel_crit]
        shell_r     = r_len[sel_crit]
        shell_pos   = self.r_pos[sel_crit]

        # getting the xyz coordinates and transforming them into spherical coordinates theta and phi:
        # axes orientation according to https://en.wikipedia.org/wiki/Spherical_coordinate_system
        if self.sphere_in_ydirection:
            shell_x = shell_pos[:,0]#2
            shell_y = shell_pos[:,1]#0
            shell_z = shell_pos[:,2]#1
            shell_phi = np.arctan2(shell_z, shell_x)
            if self.upperhalf:
                shell_theta = np.arccos(shell_y/shell_r)
            else:
                shell_theta = np.arccos(-shell_y/shell_r)
        elif self.sphere_in_xdirection:
            shell_x = shell_pos[:,0]#2
            shell_y = shell_pos[:,1]#0
            shell_z = shell_pos[:,2]#1
            shell_phi = np.arctan2(shell_y, shell_x)
            if self.upperhalf:
                shell_theta = np.arccos(shell_x/shell_r)
            else:
                shell_theta = np.arccos(-shell_x/shell_r)
        else:
            shell_x = shell_pos[:,0]#2
            shell_y = shell_pos[:,1]#0
            shell_z = shell_pos[:,2]#1
            shell_phi = np.arctan2(shell_y, shell_x)
            if self.upperhalf:
                shell_theta = np.arccos(shell_z/shell_r)
            else:
                shell_theta = np.arccos(-shell_z/shell_r)

        if volume is not None:
            if self.sphere_in_ydirection:
                if self.upperhalf:
                    vol_theta = np.arccos(self.r_pos[:,1]/r_len)
                    self.volume = volume[(r_len <= self.r) & (self.r_pos[:,1] >= self.h)]
                    self.val = value[(r_len <= self.r) & (self.r_pos[:,1] >= self.h)]
                else:
                    vol_theta = np.arccos(-self.r_pos[:,1]/r_len)
                    self.volume = volume[(r_len <= self.r) & (self.r_pos[:,1] <= -self.h)]
                    self.val = value[(r_len <= self.r) & (self.r_pos[:,1] <= -self.h)]
            elif self.sphere_in_xdirection:
                if self.upperhalf:
                    vol_theta = np.arccos(self.r_pos[:,0]/r_len)
                    self.volume = volume[(r_len <= self.r) & (self.r_pos[:,0] >= self.h)]
                    self.val = value[(r_len <= self.r) & (self.r_pos[:,0] >= self.h)]
                else:
                    vol_theta = np.arccos(-self.r_pos[:,0]/r_len)
                    self.volume = volume[(r_len <= self.r) & (self.r_pos[:,0] <= -self.h)]
                    self.val = value[(r_len <= self.r) & (self.r_pos[:,0] <= -self.h)]
            else:
                if self.upperhalf:
                    vol_theta = np.arccos(self.r_pos[:,2]/r_len)
                    self.volume = volume[(r_len <= self.r) & (self.r_pos[:,2] >= self.h)]
                    self.val = value[(r_len <= self.r) & (self.r_pos[:,2] >= self.h)]
                else:
                    vol_theta = np.arccos(-self.r_pos[:,2]/r_len)
                    self.volume = volume[(r_len <= self.r) & (self.r_pos[:,2] <= -self.h)]
                    self.val = value[(r_len <= self.r) & (self.r_pos[:,2] <= -self.h)]
            if count_jetr:
                self.counted_jetr = np.size(self.volume[self.val>1e-3])
        else:
            self.val = shell_value[shell_theta<=self.open_angle]

        self.x = shell_x[shell_theta<=self.open_angle]
        self.y = shell_y[shell_theta<=self.open_angle]
        self.z = shell_z[shell_theta<=self.open_angle]

        # self.x = shell_x[shell_y >= self.h]
        # self.y = shell_y[shell_y >= self.h]
        # self.z = shell_z[shell_y >= self.h]
        # self.val = shell_value[shell_y >= self.h]

        if self.n_top:
            np.savez('./npz/sphere_top_gridxyz.npz', data_grid_x=self.x, data_grid_y=self.y, data_grid_z=self.z)
        else:
            np.savez('./npz/sphere_bot_gridxyz.npz', data_grid_x=self.x, data_grid_y=self.y, data_grid_z=self.z)