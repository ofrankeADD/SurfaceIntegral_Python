import matplotlib.pyplot as plt
import numpy as np
from class_sphere import Sphere
from class_cone import Cone
from class_region import Region
from load_snapvalue import get_npz_data

def test_resolution_sphere_arb_field():
    # Sphere:
    h_min = 50
    h_max = 75
    dr = 5
    n_top = True
    open_angle = 30.

    Ntheta = 10
    Nphi = 10
    Nr = 10
    stepsize = 100

    path = 'p_brag_b20'

    r_len, r_pos, value_not_used = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')
    sphere = Sphere(Ntheta, Nphi, h_min, h_max, dr, n_top, r_pos, open_angle, sphere_in_ydirection=True, sphere_in_xdirection=False)

    if open_angle==30.:
        surf_int_exact = -1./960.*np.pi*sphere.r**3 * (80.*(-8.+3.*np.sqrt(3.)) + (-768. + 441.*np.sqrt(3.))*sphere.r**2)
    elif open_angle==180.:
        surf_int_exact = 4./15.*np.pi*sphere.r**3*(5. + 6.*sphere.r**2)
    else:
        print('Analytical result of surface integral for given opening angle not specified!\n')
        #return
    if not n_top:
        surf_int_exact *= -1.

    sphere.iterate_grid_resolution(tol=0.01)

    err = 1.
    tol = 0.005
    n = 1
    while err > tol:

        value = np.zeros_like(r_pos)
        value[:,0] = r_pos[:,0]**3
        value[:,1] = r_pos[:,1]
        value[:,2] = r_pos[:,2]**3

        err = sphere.iterate_surface_integral(value, surf_int_exact, stepsize, tol)
        print('Iteration: {}'.format(n))
        n += 1
    flux_sphere = sphere.flux
    print('flux_sphere: {:.2f}\n'.format(flux_sphere))

    # TODO: compare surf_int_exact here

def test_resolution_cone_arb_field():
    # Cone:
    h_min = 50
    h_max = 75
    dr = 5
    open_angle = 30.

    Ntheta = 10
    Nphi = 10
    Nr = 10
    stepsize = 100

    path = 'p_brag_b20'

    r_min = h_min*np.tan(np.radians(open_angle))
    r_max = h_max*np.tan(np.radians(open_angle))
    print(r_min, r_max)

    r_len, r_pos, value_not_used = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')

    cone = Cone(Nr, Ntheta, h_min, h_max, dr, r_pos, open_angle, cone_in_ydirection=True, cone_in_xdirection=False)
    surf_int_exact = np.pi/(30.*np.tan(cone.open_angle))*(20.*(r_max**3 - r_min**3) - 9.*(r_max**5 - r_min**5))

    cone.iterate_grid_resolution(tol=0.005)

    err = 1.
    tol = 0.005
    n = 1
    while err > tol:
        value = np.zeros_like(r_pos)
        value[:,0] = np.copy(r_pos[:,0])**3
        value[:,1] = np.copy(r_pos[:,1])
        value[:,2] = np.copy(r_pos[:,2])**3

        err = cone.iterate_surface_integral(value, surf_int_exact, stepsize, tol)
        print('Iteration: {}'.format(n))
        n += 1
    flux_cone = cone.flux
    print('flux_cone: {:.2f}\n'.format(flux_cone))

def test_resolution_sphere():
    # Sphere:
    h_min = 50
    h_max = 75
    dr = 5
    n_top = False
    open_angle = 30.

    Ntheta = 10
    Nphi = 10
    Nr = 10
    stepsize = 100

    path = 'p_brag_b20'

    r_len, r_pos, value = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')
    sphere = Sphere(Ntheta, Nphi, h_min, h_max, dr, n_top, r_pos, open_angle, sphere_in_xdirection=True)

    sphere.surface_integral(value)
    tmp_flux = sphere.flux
    sphere.iterate_grid_resolution(tol=0.01)

    err = 1.
    tol = 0.005
    n = 1
    while err > tol:
        err = sphere.iterate_surface_integral(value, tmp_flux, stepsize, tol)
        print('Iteration: {}'.format(n))
        tmp_flux = sphere.flux
        n += 1
    print('converged flux_sphere: {:.2f}\n'.format(sphere.flux))

def test_resolution_cone():
    # Cone:
    r_min = 50
    r_max = 75
    dr = 5
    open_angle = 30.

    Ntheta = 10
    Nphi = 10
    Nr = 10
    stepsize = 100

    path = 'p_brag_b20'

    r_len, r_pos, value = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')
    cone = Cone(Nr, Ntheta, r_min, r_max, dr, r_pos, open_angle, cone_in_xdirection=True)

    cone.surface_integral(value)
    tmp_flux = cone.flux
    cone.iterate_grid_resolution(tol=0.005)

    # cone.Nr += stepsize
    # cone.Ntheta += stepsize
    # cone.Nphi += stepsize
    # cone.create_grid()

    err = 1.
    tol = 0.005
    n = 1
    while err > tol:
        err = cone.iterate_surface_integral(value, tmp_flux, stepsize, tol)
        print('Iteration: {}'.format(n))
        tmp_flux = cone.flux
        n += 1
    print('converged flux_cone: {:.2f}\n'.format(cone.flux))

def test_calc_avg_surf_int_sphere():
    h_min = 50
    h_max = 75
    dr = 5
    n_top = True
    open_angle = 30.

    Ntheta = 400
    Nphi = 400
    Nr = 400
    stepsize = 100

    path = 'p_brag_b20'

    r_len, r_pos, value_not_used = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')
    sphere = Sphere(Ntheta, Nphi, h_min, h_max, dr, n_top, r_pos, open_angle, sphere_in_ydirection=True, sphere_in_xdirection=False, upperhalf=False)

    if open_angle==30.:
        surf_int_exact = -1./960.*np.pi*sphere.r**3 * (80.*(-8.+3.*np.sqrt(3.)) + (-768. + 441.*np.sqrt(3.))*sphere.r**2)
    elif open_angle==180.:
        surf_int_exact = 4./15.*np.pi*sphere.r**3*(5. + 6.*sphere.r**2)
    else:
        print('Analytical result of surface integral for given opening angle not specified!\n')
        #return
    if not n_top:
        surf_int_exact *= -1.

    print(surf_int_exact)

    value = np.zeros_like(r_pos)
    value[:,0] = r_pos[:,0]**3
    value[:,1] = r_pos[:,1]
    value[:,2] = r_pos[:,2]**3

    flux_avg = sphere.surface_integral(value, perp=False, avg=True)
    print('numerical  flux_avg:', flux_avg)
    print('analytical flux_avg:', surf_int_exact/sphere.dS_area_exact)
    #print('analytical flux_avg:', surf_int_exact/np.sum(np.sqrt(np.sum(value**2, axis=1))))
    sphere.surface_integral(value)

def test_calc_avg_surf_int_cone():
    h_min = 50
    h_max = 75
    dr = 5
    n_top = True
    open_angle = 30.

    Ntheta = 400
    Nphi = 400
    Nr = 400
    stepsize = 100

    path = 'p_brag_b20'

    r_min = h_min*np.tan(np.radians(open_angle))
    r_max = h_max*np.tan(np.radians(open_angle))

    r_len, r_pos, value_not_used = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')
    cone = Cone(Nr, Ntheta, h_min, h_max, dr, r_pos, open_angle, cone_in_ydirection=True, cone_in_xdirection=False, upperhalf=False)

    surf_int_exact = np.pi/(30.*np.tan(cone.open_angle))*(20.*(r_max**3 - r_min**3) - 9.*(r_max**5 - r_min**5))

    print(surf_int_exact)

    value = np.zeros_like(r_pos)
    value[:,0] = r_pos[:,0]**3
    value[:,1] = r_pos[:,1]
    value[:,2] = r_pos[:,2]**3

    flux_avg = cone.surface_integral(value, avg=True)
    print('numerical  flux_avg:', flux_avg)
    print('analytical flux_avg:', surf_int_exact/cone.dS_area_exact)
    #print('analytical flux_avg:', surf_int_exact/np.sum(np.sqrt(np.sum(value**2, axis=1))))
    cone.surface_integral(value)


if __name__ == "__main__":

# convergence of individual Surfaces with arbitrary vfield against known analytical solution:
    #test_resolution_sphere_arb_field()
    #test_resolution_cone_arb_field()

# convergence of individual Surfaces with velocity vector field:
    #test_resolution_sphere()
    #test_resolution_cone()

# test average surface integral flux of arbitrary vfield:
    #test_calc_avg_surf_int_sphere()
    test_calc_avg_surf_int_cone()
