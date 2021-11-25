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

