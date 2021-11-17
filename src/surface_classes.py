import numpy as np
from scipy.interpolate import NearestNDInterpolator

class Surface:
    def __init__(self):
        pass

    def iterate_surface_integral(self, value, surf_int_compare, stepsize=100, tol=0.01, perp=False):
        print('Increase resolution for object {} until flux is converged to tol={}.'.format(self.type, tol))

        self.surface_integral(value, perp)
        surf_int_approx = np.abs(self.flux)
        surf_int_cmp = np.abs(surf_int_compare)
        maxmin = np.maximum(surf_int_cmp, surf_int_approx)
        err = np.abs(surf_int_cmp - surf_int_approx)/maxmin
            
        if self.type == 'Sphere':
            print('Ntheta: {}, Nphi: {}'.format(self.Ntheta, self.Nphi))
        elif self.type == 'Cone':
            print('Nr: {}, Ntheta: {}'.format(self.Nr, self.Ntheta))
        print('step (N-1) surface integral: {:.2f}'.format(surf_int_compare))
        print('step N approx surface integral: {:.2f}'.format(self.flux))
        print('Difference: {:.2f}'.format(surf_int_compare - self.flux))
        print('Absolute error: {:.2f}%\n'.format(100.*err))

        self.Nr += stepsize
        self.Ntheta += stepsize
        self.Nphi += stepsize
        self.create_grid()

        return err

    def iterate_grid_resolution(self, stepsize=100, tol=0.01):
        # compare exact with approx surface area:
        print('Increase grid resolution for object {} until dS is converged to tol={}.'.format(self.type, tol))

        err = 1.
        n = 1
        while err > tol:
            dS = np.sqrt(self.dS_vec[:,:,0]**2 + self.dS_vec[:,:,1]**2 + self.dS_vec[:,:,2]**2)
            dS_area_approx = np.sum(dS)
            err = np.abs(self.dS_area_exact - dS_area_approx)/np.maximum(self.dS_area_exact, dS_area_approx)

            print('Iteration: {}'.format(n))
            if self.type == 'Sphere':
                print('Ntheta: {}, Nphi: {}'.format(self.Ntheta, self.Nphi))
            elif self.type == 'Cone':
                print('Nr: {}, Ntheta: {}'.format(self.Nr, self.Ntheta))
            print('Approx sphere area: {:.2f}'.format(dS_area_approx))
            print('Exact  sphere area: {:.2f}'.format(self.dS_area_exact))
            print('Difference: {:.2f}'.format(self.dS_area_exact - dS_area_approx))
            print('Absolute error: {:.2f}%'.format(100.*err))

            self.Nr += stepsize
            self.Ntheta += stepsize
            self.Nphi += stepsize
            self.create_grid()
            n += 1

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

class Cone(Surface):
    def __init__(self, Nr, Ntheta, h_min, h_max, dr, r_pos, open_angle=30., cone_in_ydirection=False, cone_in_xdirection=True, upperhalf=True):
        print('Initiating Cone object.')
        self.type = 'Cone'
        self.h_min = h_min
        self.h_max = h_max
        self.dr = dr
        self.open_angle = np.radians(open_angle)
        self.flux = 0.
        self.Nr = Nr
        self.Ntheta = Ntheta
        self.Nphi = 1
        self.r_pos = r_pos
        assert not (cone_in_ydirection==True and cone_in_xdirection==True), 'Cone object cannot be orientated in both x and y direction!'
        self.cone_in_ydirection = cone_in_ydirection
        self.cone_in_xdirection = cone_in_xdirection
        self.upperhalf = upperhalf
        self.get_dS_area_exact()
        self.get_dV_volume_exaxt()

        self.create_grid()
        Surface.__init__(self)

    def get_dS_area_exact(self):
        # Exact result from here
        # https://rechneronline.de/pi/truncated-cone.php
        r_min = self.h_min*np.tan(self.open_angle)
        r_max = self.h_max*np.tan(self.open_angle)
        R = r_max
        r = r_min
        h = self.h_max - self.h_min
        s = np.sqrt( (R-r)**2 + h**2 )
        self.dS_area_exact = (R+r)*np.pi*s

    def get_dV_volume_exaxt(self):
        r_min = self.h_min*np.tan(self.open_angle)
        r_max = self.h_max*np.tan(self.open_angle)
        R = r_max
        r = r_min
        h = self.h_max - self.h_min
        self.dV_volume_exact = h * (np.pi/3.) * (R**2 + R*r + r**2)

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

    def create_grid(self, alternative_grid=True, verbose=False):
        if verbose:
            print('Creating grid for Cone object.')

        r_min = self.h_min*np.tan(self.open_angle)
        r_max = self.h_max*np.tan(self.open_angle)

        # Resolution for the grid
        if not alternative_grid:
            dr = (r_max - r_min)/self.Nr
        else:
            dr = (r_max - r_min)/(self.Nr-1)
        dtheta = 2*np.pi/self.Ntheta

        # Grid 1D arrays
        if not alternative_grid:
            r_vec = r_min + (np.arange(self.Nr) + 1/2) * dr
            theta_vec  = (np.arange(self.Ntheta) + 1/2) * dtheta
        else:
            r_vec = np.linspace(r_min, r_max, self.Nr, endpoint=True)
            theta_vec = np.linspace(0., 2.*np.pi, self.Ntheta, endpoint=False)

        # Grid 2d arrays
        r_mat, theta_mat = np.meshgrid(r_vec, theta_vec)

        # height 2d array
        hh = r_mat/np.tan(self.open_angle)

        # Cartesian grid arrays
        self.xx = r_mat*np.cos(theta_mat)
        self.yy = r_mat*np.sin(theta_mat)
        self.zz = hh

        if not self.cone_in_xdirection:
            if self.upperhalf:
                self.zz *= -1.

        # Components of normal vector 2d array
        n_vec = np.zeros((self.Nr,self.Ntheta,3))
        n_vec[:,:,0] = -np.cos(theta_mat)/np.tan(self.open_angle)
        n_vec[:,:,1] = -np.sin(theta_mat)/np.tan(self.open_angle)
        n_vec[:,:,2] = np.ones_like(r_mat)
        self.n_vec = n_vec

        # Area element 2d arrays. \vec{dS} =  \vec{n} dx dy = \vec{n} r dr dθ
        dS_vec = np.zeros((self.Nr,self.Ntheta,3))
        dS_vec[:,:,0] = self.n_vec[:,:,0] * r_mat * dr * dtheta
        dS_vec[:,:,1] = self.n_vec[:,:,1] * r_mat * dr * dtheta
        dS_vec[:,:,2] = self.n_vec[:,:,2] * r_mat * dr * dtheta
        self.dS_vec = dS_vec

        if self.cone_in_ydirection:
            self.rotate_90deg_around_x(verbose=False)
        elif self.cone_in_xdirection:
            self.rotate_90deg_around_y(verbose=False)

        np.savez('./npz/cone_gridxxyyzz.npz', data_grid_xx=self.xx, data_grid_yy=self.yy, data_grid_zz=self.zz)

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
            print('Cone object is now orientated in y direction.')

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

        if self.cone_in_ydirection:
            if self.upperhalf:
                horiz_r = np.tan(self.open_angle) * self.r_pos[:,1]
            else:
                horiz_r = np.tan(-self.open_angle) * self.r_pos[:,1]
            d = np.sqrt(self.r_pos[:,2]**2 + self.r_pos[:,0]**2)
            h = self.r_pos[:,1]
        elif self.cone_in_xdirection:
            if self.upperhalf:
                horiz_r = np.tan(self.open_angle) * self.r_pos[:,0]
            else:
                horiz_r = np.tan(-self.open_angle) * self.r_pos[:,0]
            d = np.sqrt(self.r_pos[:,2]**2 + self.r_pos[:,1]**2)
            h = self.r_pos[:,0]
        else:
            if self.upperhalf:
                horiz_r = np.tan(self.open_angle) * self.r_pos[:,2]
            else:
                horiz_r = np.tan(-self.open_angle) * self.r_pos[:,2]
            d = np.sqrt(self.r_pos[:,1]**2 + self.r_pos[:,0]**2)
            h = self.r_pos[:,2]

        if self.upperhalf:
            sel_crit = ((h > self.h_min) & (h < self.h_max) & (horiz_r-self.dr/2. < d) & (d < horiz_r+self.dr/2.))
        else:
            sel_crit = ((h < -self.h_min) & (h > -self.h_max) & (horiz_r-self.dr/2. < d) & (d < horiz_r+self.dr/2.))
        cone_value = value[sel_crit]
        cone_r     = r_len[sel_crit]
        cone_pos   = self.r_pos[sel_crit]
        cone_x = cone_pos[:,0]
        cone_y = cone_pos[:,1]
        cone_z = cone_pos[:,2]

        if volume is not None:
            if self.upperhalf:
                self.volume = volume[(h > self.h_min) & (h < self.h_max) & (d < horiz_r)]
                self.val = value[(h > self.h_min) & (h < self.h_max) & (d < horiz_r)]
            else:
                self.volume = volume[(h < -self.h_min) & (h > -self.h_max) & (d < horiz_r)]
                self.val = value[(h < -self.h_min) & (h > -self.h_max) & (d < horiz_r)]
            if count_jetr:
                self.counted_jetr = np.size(self.volume[self.val>1e-3])
        else:
            self.val = cone_value

        self.x = cone_x
        self.y = cone_y
        self.z = cone_z
        np.savez('./npz/cone_gridxyz.npz', data_grid_x=self.x, data_grid_y=self.y, data_grid_z=self.z)

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