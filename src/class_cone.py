import numpy as np
from scipy.interpolate import NearestNDInterpolator
from class_surface import Surface

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