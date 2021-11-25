import matplotlib.pyplot as plt
import numpy as np
from class_sphere import Sphere
from class_cone import Cone
from class_region import Region
from load_snapvalue import get_npz_data

def calc_surf_int_flux_region():
    # r_min in kpc, r_max in kpc, dr in kpc, open_angle in degrees
    h_min = 50
    h_max = 75
    dr = 5
    n_top = True
    open_angle = 30.

    Ntheta = 400
    Nphi = 400
    Nr = 400

    path = 'p_brag_b20'

    r_len, r_pos, value = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue='vel')

    region = Region(Ntheta, Nphi, Nr, h_min, h_max, dr, r_pos, open_angle)
    region.get_region_surface_flux(value)
    print('total_flux: {:.2f}'.format(region.total_flux))

def calc_Q_adv_conv_sphere():
    h_min = 50
    h_max = 75
    dr = 5
    n_top = True
    open_angle = 30.

    Ntheta = 700
    Nphi = 700
    Nr = 700

    path = 'p_brag_b20'

    snapvalue = 'n'
    r_len, r_pos, n = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue=snapvalue)
    snapvalue = 'kT_erg'
    r_len, r_pos, kT = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue=snapvalue)
    snapvalue = 'vel'
    r_len, r_pos, vel = get_npz_data(use_h5py=False, rmax=100, snapnum=0, path=path, snapvalue=snapvalue)
    sphere = Sphere(Ntheta, Nphi, h_min, h_max, dr, n_top, r_pos, open_angle, sphere_in_ydirection=True, sphere_in_xdirection=False, upperhalf=True)

    flux_avg_n = sphere.surface_integral(n, avg=True)
    flux_avg_kT = sphere.surface_integral(kT, avg=True)
    flux_avg_vel = sphere.surface_integral(vel, perp=True, avg=True)
    print(flux_avg_vel.shape)

    # flux_int_n = sphere.surface_integral(n)
    # flux_int_kT = sphere.surface_integral(kT)
    # flux_int_vel = sphere.surface_integral(vel, perp=True)

    flux_del_n = n - flux_avg_n
    flux_del_kT = kT - flux_avg_kT
    flux_del_vel = vel - flux_avg_vel
    print(flux_del_vel.shape)
    del_n = np.zeros_like(vel)
    del_n[:,0] = flux_del_n
    del_n[:,1] = flux_del_n
    del_n[:,2] = flux_del_n

    del_kT = np.zeros_like(vel)
    del_kT[:,0] = flux_del_kT
    del_kT[:,1] = flux_del_kT
    del_kT[:,2] = flux_del_kT


    # sphere.get_surface_data(n)
    # # val_n = sphere.val
    # # avg_n = np.mean(val_n)
    # # del_n = val_n - avg_n
    # grid_val_n = sphere.interpolate_to_grid(perp=False)
    # grid_avg_n = np.mean(grid_val_n)
    # grid_del_n = grid_val_n - flux_avg_n
    # print(flux_avg_n)
    # print(grid_avg_n)
    # print(grid_val_n)

    # sphere.get_surface_data(kT)
    # # val_kT = sphere.val
    # # avg_kT = np.mean(val_kT)
    # # del_kT = val_kT - avg_kT
    # grid_val_kT = sphere.interpolate_to_grid(perp=False)
    # grid_avg_kT = np.mean(grid_val_kT)
    # grid_del_kT = grid_val_kT - flux_avg_kT
    # print(flux_avg_kT)
    # print(grid_avg_kT)
    # print(grid_val_kT)

    # sphere.get_surface_data(vel)
    # # val_vel = sphere.val
    # # avg_vel = np.mean(val_vel)
    # # del_vel = val_vel - avg_vel
    # grid_val_vel = sphere.interpolate_to_grid(perp=True)
    # grid_avg_vel = np.mean(grid_val_vel)
    # grid_del_vel = grid_val_vel - flux_avg_vel
    # print(flux_avg_vel)
    # print(grid_avg_vel)
    # print(grid_val_vel)


    ##Q_adv  = 1.5 * (grid_avg_n*grid_avg_kT*grid_avg_vel + grid_avg_kT*np.mean(grid_del_n*grid_del_vel))
    Q_adv  = 1.5 * (flux_avg_n*flux_avg_kT*flux_avg_vel + flux_avg_kT*sphere.surface_integral(del_n*flux_del_vel, perp=True, avg=True))
    ##Q_adv  = 1.5 * (flux_avg_n*flux_avg_kT*flux_avg_vel + flux_avg_kT*sphere.calculate_avg_surf_int(flux_del_n)*sphere.calculate_avg_surf_int(flux_del_vel, perp=True))
    #Q_adv  = 1.5 * (flux_avg_n*flux_avg_kT*flux_avg_vel + flux_avg_kT*np.mean(grid_del_n*grid_del_vel))
    
    ##Q_conv = 1.5 * (grid_avg_n*np.mean(grid_del_kT*grid_del_vel) + grid_avg_vel*np.mean(grid_del_n*grid_del_kT) + np.mean(grid_del_n*grid_del_kT*grid_del_vel))
    Q_conv = 1.5 * (flux_avg_n*sphere.surface_integral(del_kT*flux_del_vel, perp=True, avg=True) + flux_avg_vel*sphere.surface_integral(flux_del_n*flux_del_kT, perp=False, avg=True) + sphere.surface_integral(del_n*del_kT*flux_del_vel, perp=True, avg=True))
    ##Q_conv = 1.5 * (flux_avg_n*sphere.calculate_avg_surf_int(flux_del_kT)*sphere.calculate_avg_surf_int(flux_del_vel, perp=True) + flux_avg_vel*sphere.calculate_avg_surf_int(flux_del_n)*sphere.calculate_avg_surf_int(flux_del_kT) + sphere.calculate_avg_surf_int(flux_del_n)*sphere.calculate_avg_surf_int(flux_del_kT)*sphere.calculate_avg_surf_int(flux_del_vel, perp=True))
    #Q_conv = 1.5 * (flux_avg_n*np.mean(grid_del_kT*grid_del_vel) + flux_avg_vel*np.mean(grid_del_n*grid_del_kT) + np.mean(grid_del_n*grid_del_kT*grid_del_vel))

    flux_int_Q_adv = sphere.surface_integral(np.ones_like(n)*Q_adv, perp=False, avg=False)
    print(Q_adv*sphere.area)
    flux_int_Q_conv = sphere.surface_integral(np.ones_like(n)*Q_conv, perp=False, avg=False)
    print(Q_conv*sphere.area)

def calc_L_adv_conv_jet_region():
    import gc
    import subprocess
    import matplotlib as mpl
    mpl.use('Agg')

    kpc_in_cm = 3.085677581e21
    UnitLength_in_cm = 3.085678e24 #  1.0 Mpc

    dr = 3
    dbin = 8
    rmax = 100
    open_angle = 30.
    Ntheta = 400
    Nphi = 400
    Nr = 400

    snapmax = 85
    if snapmax == 85:
        makemovie = True
    else:
        makemovie = False

    plot_grid = False
    count_jetr = True

    binmin = dbin+dr/2.
    binmax = rmax/(2. - np.cos(np.radians(open_angle)))
    bins = np.arange(binmin,binmax,dbin)
    print(bins)
    path = 'p_brag_b20'
    filename = 'bins_pn_jet_qpqm_H_adv_conv_bins{}_'.format(dbin)+path+'_'
    print(filename)


    clear_npz = False
    if clear_npz:
        subprocess.call('rm ./npz/'+filename+'*.npz', shell=True)

    try:
        npz_mean = np.load('./npz/'+filename+'mean.npz')
        mean_qp_exact = npz_mean['data_mean_qp_exact']
        mean_qp_approx = npz_mean['data_mean_qp_approx']
        mean_qm_exact = npz_mean['data_mean_qm_exact']
        mean_qm_approx = npz_mean['data_mean_qm_approx']
        mean_H_adv = npz_mean['data_mean_H_adv']
        mean_H_conv = npz_mean['data_mean_H_conv']
        mean_count = npz_mean['data_count']
    except:
        mean_qp_exact = np.zeros((len(bins)))
        mean_qp_approx = np.zeros((len(bins)))
        mean_qm_exact = np.zeros((len(bins)))
        mean_qm_approx = np.zeros((len(bins)))
        mean_H_adv = np.zeros((len(bins)))
        mean_H_conv = np.zeros((len(bins)))
        mean_count = 0

        # loop over snapshots:
        for snapn in range(0, snapmax, 1):
            if (snapn%2)==0 or (snapn%2)!=0:
                print('SNAPN:', snapn)

                snapvalue = 'n'
                r_len, r_pos, time_myr, n = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'kT_erg'
                r_len, r_pos, time_myr, kT = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'vel_cgs'
                r_len, r_pos, time_myr, vel = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'vol'
                r_len, r_pos, time_myr, volume = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'qp_data'
                r_len, r_pos, time_myr, qp = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'qm_theo'
                r_len, r_pos, time_myr, qm = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                if count_jetr:
                    snapvalue = 'jetr'
                    r_len, r_pos, time_myr, jetr = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)

                bins_h_dr = []
                bins_H_adv = []
                bins_H_adv_pos = []
                bins_H_adv_neg = []
                bins_H_conv = []
                bins_H_conv_pos = []
                bins_H_conv_neg = []
                bins_vol_exact = []
                bins_vol_approx = []
                bins_qp = []
                bins_qm = []
                tmp_L_adv_top = np.zeros_like(bins)
                tmp_L_conv_top = np.zeros_like(bins)

                try:
                    npz = np.load('./npz/'+filename+'{:02d}.npz'.format(snapn))
                    bins_h_dr = npz['data_bin']
                    bins_H_adv = npz['data_H_adv']
                    bins_H_adv_pos = npz['data_H_adv_pos']
                    bins_H_adv_neg = npz['data_H_adv_neg']
                    bins_H_conv = npz['data_H_conv']
                    bins_H_conv_pos = npz['data_H_conv_pos']
                    bins_H_conv_neg = npz['data_H_conv_neg']
                    bins_vol_exact = npz['data_vol_exact']
                    bins_vol_approx = npz['data_vol_approx']
                    bins_qp = npz['data_qp']
                    bins_qm = npz['data_qm']
                except:

                    # loop over cone segments:
                    for i,b in enumerate(bins):
                        h_min = (b-dbin)
                        h_max = b

                        region = Region(Ntheta, Nphi, Nr, h_min*kpc_in_cm, h_max*kpc_in_cm, dr*kpc_in_cm, r_pos*kpc_in_cm, open_angle, region_in_ydirection=False, region_in_xdirection=True, upperhalf=True)
                        print('bin:', b)
                        print('h_min:', h_min)
                        print('h_max:', h_max)

                        L_adv_top, L_adv_bot, L_adv_con, L_conv_top, L_conv_bot, L_conv_con = region.calc_Q_adv_conv_region(n, kT, vel)
                        if i > 0:
                            assert L_adv_bot == -tmp_L_adv_top[i-1], 'Spherical bottom cap should have the same absolute adv flux as the top cap from the previuos segment bin!'
                            assert L_conv_bot == -tmp_L_conv_top[i-1], 'Spherical bottom cap should have the same absolute conv flux as the top cap from the previuos segment bin!'
                            #print(L_adv_bot, tmp_L_adv_top[i-1])
                            #print(L_conv_bot, tmp_L_conv_top[i-1])
                        tmp_L_adv_top[i] = L_adv_top
                        tmp_L_conv_top[i] = L_conv_top

                        print('L_adv_top:', L_adv_top)
                        print('L_adv_bot:', L_adv_bot)
                        print('L_adv_con:', L_adv_con)
                        print('L_conv_top:', L_conv_top)
                        print('L_conv_bot:', L_conv_bot)
                        print('L_conv_con:', L_conv_con)

                        L_adv = L_adv_top + L_adv_bot + L_adv_con
                        L_conv = L_conv_top + L_conv_bot + L_conv_con

                        if plot_grid and snapn == 0:
                            region.plot_xxyyzz_grid(snapn, b, kpc_in_cm)
                            region.plot_xyz_grid(snapn, b, kpc_in_cm)

                        dV_exact_top, dV_exact_bot, dV_exact_con = region.get_region_dV_volume_exact()
                        if count_jetr:
                            dV_approx_top, dV_approx_bot, dV_approx_con, dV_jetr_top, dV_jetr_bot, dV_jetr_con = region.get_region_dV_volume_approx(jetr, volume*UnitLength_in_cm**3, count_jetr=True)
                        dV_approx_top, dV_approx_bot, dV_approx_con, dV_qp_top, dV_qp_bot, dV_qp_con = region.get_region_dV_volume_approx(qp, volume*UnitLength_in_cm**3)
                        dV_approx_top, dV_approx_bot, dV_approx_con, dV_qm_top, dV_qm_bot, dV_qm_con = region.get_region_dV_volume_approx(qm, volume*UnitLength_in_cm**3)
                        dV_region_exact = dV_exact_top - dV_exact_bot + dV_exact_con
                        dV_region_approx = dV_approx_top - dV_approx_bot + dV_approx_con
                        print(dV_region_exact)
                        print(dV_region_approx)

                        H_adv = L_adv/dV_region_approx
                        #H_adv_con = (L_adv_con)/(dV_approx_con)
                        H_conv = L_conv/dV_region_approx
                        #H_conv_con = L_conv_con/dV_approx_con
                        print('H_adv:', H_adv)
                        #print('H_adv_con:', H_adv_con)
                        print('H_conv:', H_conv)
                        #print('H_conv_con:', H_conv_con)

                        bins_vol_exact.append([dV_exact_top, dV_exact_bot, dV_exact_con])
                        bins_vol_approx.append([dV_approx_top, dV_approx_bot, dV_approx_con])
                        bins_h_dr.append((h_min+h_max)/2-dr)

                        bins_H_adv.append([L_adv_top, L_adv_bot, L_adv_con])
                        bins_H_conv.append([L_conv_top, L_conv_bot, L_conv_con])

                        bins_qp.append([dV_qp_top, dV_qp_bot, dV_qp_con])
                        bins_qm.append([dV_qm_top, dV_qm_bot, dV_qm_con])

                        if H_adv >= 0.:
                            bins_H_adv_pos.append([(h_min+h_max)/2-dr, H_adv])
                        else:
                            bins_H_adv_neg.append([(h_min+h_max)/2-dr, -H_adv])
                        if H_conv >= 0.:
                            bins_H_conv_pos.append([(h_min+h_max)/2-dr, H_conv])
                        else:
                            bins_H_conv_neg.append([(h_min+h_max)/2-dr, -H_conv])

                        del region
                        gc.collect()
                
                    np.savez('./npz/'+filename+'{:02d}.npz'.format(snapn), data_bin=bins_h_dr, data_H_adv=bins_H_adv, data_H_conv=bins_H_conv, data_vol_exact=bins_vol_exact, data_vol_approx=bins_vol_approx, data_H_adv_pos=bins_H_adv_pos, data_H_adv_neg=bins_H_adv_neg, data_H_conv_pos=bins_H_conv_pos, data_H_conv_neg=bins_H_conv_neg, data_qp=bins_qp, data_qm=bins_qm)
                
                adv = (np.array(bins_H_adv)[:,0]+np.array(bins_H_adv)[:,1]+np.array(bins_H_adv)[:,2])
                conv = (np.array(bins_H_conv)[:,0]+np.array(bins_H_conv)[:,1]+np.array(bins_H_conv)[:,2])
                qp = (np.array(bins_qp)[:,0]-np.array(bins_qp)[:,1]+np.array(bins_qp)[:,2])
                qm = (np.array(bins_qm)[:,0]-np.array(bins_qm)[:,1]+np.array(bins_qm)[:,2])
                dV_approx = (np.array(bins_vol_approx)[:,0]-np.array(bins_vol_approx)[:,1]+np.array(bins_vol_approx)[:,2])
                dV_exact = (np.array(bins_vol_exact)[:,0]-np.array(bins_vol_exact)[:,1]+np.array(bins_vol_exact)[:,2])

                adv_approx = adv / dV_approx
                conv_approx = conv / dV_approx
                mean_H_adv += adv_approx
                mean_H_conv += conv_approx
                
                qp_exact = qp / dV_exact
                qp_approx = qp / dV_approx
                qm_exact = qm / dV_exact
                qm_approx = qm / dV_approx
                mean_qp_exact += qp_exact
                mean_qp_approx += qp_approx
                mean_qm_exact += qm_exact
                mean_qm_approx += qm_approx
                mean_count += 1

                fig = plt.figure(num=1)
                plt.clf()
                fig, ax = plt.subplots(num=1)

                ax.plot(bins_h_dr, np.abs(qm_approx), color='C0', ls='--', label=r'$Q^-$')
                ax.plot(bins_h_dr, np.abs(qp_approx), color='C1', ls='-', label=r'$Q^+_{\mathrm{lim}}$')

                if len(bins_H_adv_pos) > 0:
                    ax.scatter(list(np.array(bins_H_adv_pos)[:,0]),  list(np.array(bins_H_adv_pos)[:,1]), color='r', marker='+', label=r'$H^+_{adv}$')
                if len(bins_H_adv_neg) > 0:
                    ax.scatter(list(np.array(bins_H_adv_neg)[:,0]),  list(np.array(bins_H_adv_neg)[:,1]), color='b', marker='+', label=r'$H^-_{adv}$')
                if len(bins_H_conv_pos) > 0:
                    ax.scatter(list(np.array(bins_H_conv_pos)[:,0]), list(np.array(bins_H_conv_pos)[:,1]), facecolors='none', edgecolors='r', label=r'$H^+_{conv}$')
                if len(bins_H_conv_neg) > 0:
                    ax.scatter(list(np.array(bins_H_conv_neg)[:,0]), list(np.array(bins_H_conv_neg)[:,1]), facecolors='none', edgecolors='b', label=r'$H^-_{conv}$')
                
                ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
                ax.set_ylabel(r'$H_{\mathrm{jet}},C_{\mathrm{jet}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
                ax.set_yscale('log')
                ax.set_xlim([0,binmax])
                ax.set_ylim([1e-28,1e-22])
                ax.legend()
                ax.set_title(r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
                fig.tight_layout()
                plt.savefig('./images/'+filename+'{:02d}.png'.format(snapn), dpi=300)
                plt.close('all')

        np.savez('./npz/'+filename+'mean.npz', data_mean_qp_exact=mean_qp_exact, data_mean_qp_approx=mean_qp_approx, data_mean_qm_exact=mean_qm_exact, data_mean_qm_approx=mean_qm_approx, data_mean_H_adv=mean_H_adv, data_mean_H_conv=mean_H_conv, data_count=mean_count)

    fig = plt.figure(num=1)
    plt.clf()
    fig, ax = plt.subplots(num=1)

    npz = np.load('./npz/'+filename+'00.npz')
    bins_h_dr = npz['data_bin']
    #ax.plot(bins_h_dr, mean_qm_exact/mean_count, color='k', ls='--', label=r'$<Q^->$')
    ax.plot(bins_h_dr, np.abs(mean_qm_approx)/mean_count, color='C0', ls='--', label=r'$<Q^->$')
    #ax.plot(bins_h_dr, mean_qp_exact/mean_count, color='k', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(bins_h_dr, np.abs(mean_qp_approx)/mean_count, color='C1', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(bins_h_dr, np.abs(mean_H_adv)/mean_count, color='C2', ls='dotted', label=r'$<H_{\mathrm{adv}}>$')
    ax.plot(bins_h_dr, np.abs(mean_H_conv)/mean_count, color='C3', ls='dashdot', label=r'$<H_{\mathrm{conv}}>$')
    
    ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
    ax.set_ylabel(r'$H_{\mathrm{jet}},C_{\mathrm{jet}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
    ax.set_yscale('log')
    ax.set_xlim([0,binmax])
    ax.set_ylim([1e-28,1e-22])
    ax.legend()
    snapvalue = 'n'
    r_len, r_pos, time_myr, n = get_npz_data(use_h5py=False, rmax=rmax, snapnum=mean_count, path=path, snapvalue=snapvalue)  
    ax.set_title('mean over time '+r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
    fig.tight_layout()
    plt.savefig('./images/'+filename+'mean.png', dpi=300)
    plt.close('all')

    if makemovie:
        subprocess.call('ffmpeg -y -framerate 10 -i ./images/'+filename+'%2d.png -pix_fmt yuv420p -vf scale=1280:-2 ./movies/'+filename+'movie.mp4', shell=True)
        subprocess.call('ffmpeg -y -i ./movies/'+filename+'movie.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 ./movies/'+filename+'movie.gif', shell=True)

def calc_L_adv_conv_shell_region():
    import gc
    import subprocess
    import matplotlib as mpl
    mpl.use('Agg')

    kpc_in_cm = 3.085677581e21
    UnitLength_in_cm = 3.085678e24 #  1.0 Mpc

    dr = 3
    dbin = 8
    rmax = 100
    open_angle = 90.
    cone_angle = 30.
    Ntheta = 400
    Nphi = 400
    Nr = 400

    snapmax = 85
    if snapmax == 85:
        makemovie = True
    else:
        makemovie = False

    plot_grid = False
    count_jetr = True

    binmin = dbin+dr/2.
    binmax = rmax/(2. - np.cos(np.radians(cone_angle)))
    bins = np.arange(binmin,binmax,dbin)
    print(bins)
    #bins = np.linspace(binmin,binmax,dbin, endpoint=True)
    #print(bins)
    path = 'p_brag_b20'
    filename = 'bins_pn_all_qpqm_H_adv_conv_bins{}_'.format(dbin)+path+'_'


    clear_npz = False
    if clear_npz:
        subprocess.call('rm ./npz/'+filename+'*.npz', shell=True)

    try:
        npz_mean = np.load('./npz/'+filename+'mean.npz')
        mean_qp_exact = npz_mean['data_mean_qp_exact']
        mean_qp_approx = npz_mean['data_mean_qp_approx']
        mean_qm_exact = npz_mean['data_mean_qm_exact']
        mean_qm_approx = npz_mean['data_mean_qm_approx']
        mean_H_adv = npz_mean['data_mean_H_adv']
        mean_H_conv = npz_mean['data_mean_H_conv']
        mean_count = npz_mean['data_count']
    except:
        mean_qp_exact = np.zeros((len(bins)))
        mean_qp_approx = np.zeros((len(bins)))
        mean_qm_exact = np.zeros((len(bins)))
        mean_qm_approx = np.zeros((len(bins)))
        mean_H_adv = np.zeros((len(bins)))
        mean_H_conv = np.zeros((len(bins)))
        mean_count = 0

        # loop over snapshots:
        for snapn in range(0, snapmax, 1):
            if (snapn%2)==0 or (snapn%2)!=0:
                print('SNAPN:', snapn)

                snapvalue = 'n'
                r_len, r_pos, time_myr, n = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'kT_erg'
                r_len, r_pos, time_myr, kT = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'vel_cgs'
                r_len, r_pos, time_myr, vel = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'vol'
                r_len, r_pos, time_myr, volume = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'qp_data'
                r_len, r_pos, time_myr, qp = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                snapvalue = 'qm_theo'
                r_len, r_pos, time_myr, qm = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)
                if count_jetr:
                    snapvalue = 'jetr'
                    r_len, r_pos, time_myr, jetr = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue=snapvalue)

                bins_h_dr = []
                bins_H_adv = []
                bins_H_adv_pos = []
                bins_H_adv_neg = []
                bins_H_conv = []
                bins_H_conv_pos = []
                bins_H_conv_neg = []
                bins_vol_exact = []
                bins_vol_approx = []
                bins_qp = []
                bins_qm = []
                tmp_L_adv_top = np.zeros_like(bins)
                tmp_L_conv_top = np.zeros_like(bins)

                try:
                    npz = np.load('./npz/'+filename+'{:02d}.npz'.format(snapn))
                    bins_h_dr = npz['data_bin']
                    bins_H_adv = npz['data_H_adv']
                    bins_H_adv_pos = npz['data_H_adv_pos']
                    bins_H_adv_neg = npz['data_H_adv_neg']
                    bins_H_conv = npz['data_H_conv']
                    bins_H_conv_pos = npz['data_H_conv_pos']
                    bins_H_conv_neg = npz['data_H_conv_neg']
                    bins_vol_exact = npz['data_vol_exact']
                    bins_vol_approx = npz['data_vol_approx']
                    bins_qp = npz['data_qp']
                    bins_qm = npz['data_qm']
                except:

                    # loop over cone segments:
                    for i,b in enumerate(bins):
                        h_min = (b-dbin)
                        h_max = b

                        region = Region(Ntheta, Nphi, Nr, h_min*kpc_in_cm, h_max*kpc_in_cm, dr*kpc_in_cm, r_pos*kpc_in_cm, open_angle, region_in_ydirection=False, region_in_xdirection=True, upperhalf=True, cone_angle=cone_angle)
                        print('bin:', b)
                        print('h_min:', h_min)
                        print('h_max:', h_max)

                        L_adv_top, L_adv_bot, L_adv_con, L_conv_top, L_conv_bot, L_conv_con = region.calc_Q_adv_conv_region(n, kT, vel)
                        if i > 0:
                            assert L_adv_bot == -tmp_L_adv_top[i-1], 'Spherical bottom cap should have the same absolute adv flux as the top cap from the previuos segment bin!'
                            assert L_conv_bot == -tmp_L_conv_top[i-1], 'Spherical bottom cap should have the same absolute conv flux as the top cap from the previuos segment bin!'
                            #print(L_adv_bot, tmp_L_adv_top[i-1])
                            #print(L_conv_bot, tmp_L_conv_top[i-1])
                        tmp_L_adv_top[i] = L_adv_top
                        tmp_L_conv_top[i] = L_conv_top

                        print('L_adv_top:', L_adv_top)
                        print('L_adv_bot:', L_adv_bot)
                        print('L_conv_top:', L_conv_top)
                        print('L_conv_bot:', L_conv_bot)

                        L_adv = L_adv_top + L_adv_bot
                        L_conv = L_conv_top + L_conv_bot

                        if plot_grid and snapn == 0:
                            region.plot_xxyyzz_grid(snapn, b, kpc_in_cm)
                            region.plot_xyz_grid(snapn, b, kpc_in_cm)

                        dV_exact_top, dV_exact_bot, dV_exact_con = region.get_region_dV_volume_exact()
                        if count_jetr:
                            dV_approx_top, dV_approx_bot, dV_approx_con, dV_jetr_top, dV_jetr_bot, dV_jetr_con = region.get_region_dV_volume_approx(jetr, volume*UnitLength_in_cm**3, count_jetr=True)
                        dV_approx_top, dV_approx_bot, dV_approx_con, dV_qp_top, dV_qp_bot, dV_qp_con = region.get_region_dV_volume_approx(qp, volume*UnitLength_in_cm**3)
                        dV_approx_top, dV_approx_bot, dV_approx_con, dV_qm_top, dV_qm_bot, dV_qm_con = region.get_region_dV_volume_approx(qm, volume*UnitLength_in_cm**3)
                        dV_region_exact = dV_exact_top - dV_exact_bot
                        dV_region_approx = dV_approx_top - dV_approx_bot
                        print(dV_region_exact)
                        print(dV_region_approx)

                        H_adv = L_adv/dV_region_approx
                        H_adv_shell = L_adv_top/dV_approx_top + L_adv_bot/dV_approx_bot
                        H_conv = L_conv/dV_region_approx
                        print('H_adv:', H_adv)
                        print('H_conv:', H_conv)

                        bins_vol_exact.append([dV_exact_top, dV_exact_bot])
                        bins_vol_approx.append([dV_approx_top, dV_approx_bot])
                        bins_h_dr.append((h_min+h_max)/2-dr)

                        bins_H_adv.append([L_adv_top, L_adv_bot])
                        bins_H_conv.append([L_conv_top, L_conv_bot])

                        bins_qp.append([dV_qp_top, dV_qp_bot])
                        bins_qm.append([dV_qm_top, dV_qm_bot])

                        if H_adv >= 0.:
                            bins_H_adv_pos.append([(h_min+h_max)/2-dr, H_adv])
                        else:
                            bins_H_adv_neg.append([(h_min+h_max)/2-dr, -H_adv])
                        if H_conv >= 0.:
                            bins_H_conv_pos.append([(h_min+h_max)/2-dr, H_conv])
                        else:
                            bins_H_conv_neg.append([(h_min+h_max)/2-dr, -H_conv])

                        del region
                        gc.collect()
                
                    np.savez('./npz/'+filename+'{:02d}.npz'.format(snapn), data_bin=bins_h_dr, data_H_adv=bins_H_adv, data_H_conv=bins_H_conv, data_vol_exact=bins_vol_exact, data_vol_approx=bins_vol_approx, data_H_adv_pos=bins_H_adv_pos, data_H_adv_neg=bins_H_adv_neg, data_H_conv_pos=bins_H_conv_pos, data_H_conv_neg=bins_H_conv_neg, data_qp=bins_qp, data_qm=bins_qm)
                
                adv = (np.array(bins_H_adv)[:,0]+np.array(bins_H_adv)[:,1])
                conv = (np.array(bins_H_conv)[:,0]+np.array(bins_H_conv)[:,1])
                qp = (np.array(bins_qp)[:,0]-np.array(bins_qp)[:,1])
                qm = (np.array(bins_qm)[:,0]-np.array(bins_qm)[:,1])
                dV_approx = (np.array(bins_vol_approx)[:,0]-np.array(bins_vol_approx)[:,1])
                dV_exact = (np.array(bins_vol_exact)[:,0]-np.array(bins_vol_exact)[:,1])

                adv_approx = adv / dV_approx
                conv_approx = conv / dV_approx
                mean_H_adv += adv_approx
                mean_H_conv += conv_approx

                qp_exact = qp / dV_exact
                qp_approx = qp / dV_approx
                qm_exact = qm / dV_exact
                qm_approx = qm / dV_approx
                mean_qp_exact += qp_exact
                mean_qp_approx += qp_approx
                mean_qm_exact += qm_exact
                mean_qm_approx += qm_approx
                mean_count += 1

                fig = plt.figure(num=1)
                plt.clf()
                fig, ax = plt.subplots(num=1)

                ax.plot(bins_h_dr, np.abs(qm_approx), color='C0', ls='--', label=r'$Q^-$')
                ax.plot(bins_h_dr, np.abs(qp_approx), color='C1', ls='-', label=r'$Q^+_{\mathrm{lim}}$')

                if len(bins_H_adv_pos) > 0:
                    ax.scatter(list(np.array(bins_H_adv_pos)[:,0]),  list(np.array(bins_H_adv_pos)[:,1]), color='r', marker='+', label=r'$H^+_{adv}$')
                if len(bins_H_adv_neg) > 0:
                    ax.scatter(list(np.array(bins_H_adv_neg)[:,0]),  list(np.array(bins_H_adv_neg)[:,1]), color='b', marker='+', label=r'$H^-_{adv}$')
                if len(bins_H_conv_pos) > 0:
                    ax.scatter(list(np.array(bins_H_conv_pos)[:,0]), list(np.array(bins_H_conv_pos)[:,1]), facecolors='none', edgecolors='r', label=r'$H^+_{conv}$')
                if len(bins_H_conv_neg) > 0:
                    ax.scatter(list(np.array(bins_H_conv_neg)[:,0]), list(np.array(bins_H_conv_neg)[:,1]), facecolors='none', edgecolors='b', label=r'$H^-_{conv}$')
                
                ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
                ax.set_ylabel(r'$H_{\mathrm{all}},C_{\mathrm{all}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
                ax.set_yscale('log')
                ax.set_xlim([0,binmax])
                ax.set_ylim([1e-28,1e-22])
                ax.legend()
                ax.set_title(r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
                fig.tight_layout()
                plt.savefig('./images/'+filename+'{:02d}.png'.format(snapn), dpi=300)
                plt.close('all')

        np.savez('./npz/'+filename+'mean.npz', data_mean_qp_exact=mean_qp_exact, data_mean_qp_approx=mean_qp_approx, data_mean_qm_exact=mean_qm_exact, data_mean_qm_approx=mean_qm_approx, data_mean_H_adv=mean_H_adv, data_mean_H_conv=mean_H_conv, data_count=mean_count)

    fig = plt.figure(num=1)
    plt.clf()
    fig, ax = plt.subplots(num=1)

    npz = np.load('./npz/'+filename+'00.npz')
    bins_h_dr = npz['data_bin']
    #ax.plot(bins_h_dr, mean_qm_exact/mean_count, color='k', ls='--', label=r'$<Q^->$')
    ax.plot(bins_h_dr, np.abs(mean_qm_approx)/mean_count, color='C0', ls='--', label=r'$<Q^->$')
    #ax.plot(bins_h_dr, mean_qp_exact/mean_count, color='k', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(bins_h_dr, np.abs(mean_qp_approx)/mean_count, color='C1', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(bins_h_dr, np.abs(mean_H_adv)/mean_count, color='C2', ls='dotted', label=r'$<H_{\mathrm{adv}}>$')
    ax.plot(bins_h_dr, np.abs(mean_H_conv)/mean_count, color='C3', ls='dashdot', label=r'$<H_{\mathrm{conv}}>$')
    
    ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
    ax.set_ylabel(r'$H_{\mathrm{all}},C_{\mathrm{all}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
    ax.set_yscale('log')
    ax.set_xlim([0,binmax])
    ax.set_ylim([1e-28,1e-22])
    ax.legend()
    snapvalue = 'n'
    r_len, r_pos, time_myr, n = get_npz_data(use_h5py=False, rmax=rmax, snapnum=mean_count, path=path, snapvalue=snapvalue)      
    ax.set_title('mean over time '+r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
    fig.tight_layout()
    plt.savefig('./images/'+filename+'mean.png', dpi=300)
    plt.close('all')

    if makemovie:
        subprocess.call('ffmpeg -y -framerate 10 -i ./images/'+filename+'%2d.png -pix_fmt yuv420p -vf scale=1280:-2 ./movies/'+filename+'movie.mp4', shell=True)
        subprocess.call('ffmpeg -y -i ./movies/'+filename+'movie.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 ./movies/'+filename+'movie.gif', shell=True)

def calc_L_adv_conv_ambient_region():
    import gc
    import subprocess
    import matplotlib as mpl
    mpl.use('Agg')

    dr = 3
    dbin = 8
    rmax = 100
    open_angle = 30.

    snapmax = 85
    if snapmax == 85:
        makemovie = True
    else:
        makemovie = False

    binmin = dbin+dr/2.
    binmax = rmax/(2. - np.cos(np.radians(open_angle)))
    bins = np.arange(binmin,binmax,dbin)
    print(bins)

    path = 'p_brag_b20'
    filename_all = 'bins_pn_all_qpqm_H_adv_conv_bins{}_'.format(dbin)+path+'_'
    filename_jet = 'bins_pn_jet_qpqm_H_adv_conv_bins{}_'.format(dbin)+path+'_'
    filename_amb = 'bins_pn_amb_qpqm_H_adv_conv_bins{}_'.format(dbin)+path+'_'


    clear_npz = True
    if clear_npz:
        ##subprocess.call('rm ./npz/'+filename_all+'mean.npz', shell=True)
        ##subprocess.call('rm ./npz/'+filename_jet+'mean.npz', shell=True)
        subprocess.call('rm ./npz/'+filename_amb+'mean.npz', shell=True)
        subprocess.call('rm ./npz/'+filename_amb+'*.npz', shell=True)

    try:
        npz_mean_all = np.load('./npz/'+filename_all+'mean.npz')
        #all_mean_qp_exact = npz_mean_all['data_mean_qp_exact']
        all_mean_qp_approx = npz_mean_all['data_mean_qp_approx']
        #all_mean_qm_exact = npz_mean_all['data_mean_qm_exact']
        all_mean_qm_approx = npz_mean_all['data_mean_qm_approx']
        all_mean_H_adv = npz_mean_all['data_mean_H_adv']
        all_mean_H_conv = npz_mean_all['data_mean_H_conv']
        all_mean_count = npz_mean_all['data_count']

        npz_mean_jet = np.load('./npz/'+filename_jet+'mean.npz')
        #jet_mean_qp_exact = npz_mean_jet['data_mean_qp_exact']
        jet_mean_qp_approx = npz_mean_jet['data_mean_qp_approx']
        #jet_mean_qm_exact = npz_mean_jet['data_mean_qm_exact']
        jet_mean_qm_approx = npz_mean_jet['data_mean_qm_approx']
        jet_mean_H_adv = npz_mean_jet['data_mean_H_adv']
        jet_mean_H_conv = npz_mean_jet['data_mean_H_conv']
        jet_mean_count = npz_mean_jet['data_count']

        npz_mean_amb = np.load('./npz/'+filename_amb+'mean.npz')
        #amb_mean_qp_exact = npz_mean_amb['data_mean_qp_exact']
        amb_mean_qp_approx = npz_mean_amb['data_mean_qp_approx']
        #amb_mean_qm_exact = npz_mean_amb['data_mean_qm_exact']
        amb_mean_qm_approx = npz_mean_amb['data_mean_qm_approx']
        amb_mean_H_adv = npz_mean_amb['data_mean_H_adv']
        amb_mean_H_conv = npz_mean_amb['data_mean_H_conv']
        amb_mean_count = npz_mean_amb['data_count']
    except:

        all_mean_qp_approx = np.zeros((len(bins)))
        all_mean_qm_approx = np.zeros((len(bins)))
        all_mean_H_adv = np.zeros((len(bins)))
        all_mean_H_conv = np.zeros((len(bins)))
        all_mean_count = 0

        jet_mean_qp_approx = np.zeros((len(bins)))
        jet_mean_qm_approx = np.zeros((len(bins)))
        jet_mean_H_adv = np.zeros((len(bins)))
        jet_mean_H_conv = np.zeros((len(bins)))
        jet_mean_count = 0

        amb_mean_qp_approx = np.zeros((len(bins)))
        amb_mean_qm_approx = np.zeros((len(bins)))
        amb_mean_H_adv = np.zeros((len(bins)))
        amb_mean_H_conv = np.zeros((len(bins)))
        amb_mean_count = 0

        for snapn in range(0, snapmax, 1):
            if snapn != 87:
                print('SNAPN:', snapn)

                # just to get the time_myr:
                r_len, r_pos, time_myr, n = get_npz_data(use_h5py=False, rmax=rmax, snapnum=snapn, path=path, snapvalue='n')

                try:
                    npz_amb = np.load('./npz/'+filename_amb+'{:02d}.npz'.format(snapn))
                    amb_bins_h_dr = npz_amb['data_bin']
                    amb_adv_approx = npz_amb['data_H_adv']
                    amb_conv_approx = npz_amb['data_H_conv']
                    amb_qp_approx = npz_amb['data_H_qp']
                    amb_qm_approx = npz_amb['data_H_qm']
                    amb_bins_H_adv_pos = npz_amb['data_H_adv_pos']
                    amb_bins_H_adv_neg = npz_amb['data_H_adv_neg']
                    amb_bins_H_conv_pos = npz_amb['data_H_conv_pos']
                    amb_bins_H_conv_neg = npz_amb['data_H_conv_neg']
                except:
                    amb_bins_H_adv_pos = []
                    amb_bins_H_adv_neg = []
                    amb_bins_H_conv_pos = []
                    amb_bins_H_conv_neg = []

                    npz_all = np.load('./npz/'+filename_all+'{:02d}.npz'.format(snapn))
                    all_bins_h_dr = npz_all['data_bin']
                    all_bins_H_adv = npz_all['data_H_adv']
                    all_bins_H_conv = npz_all['data_H_conv']
                    all_bins_vol_approx = npz_all['data_vol_approx']
                    all_bins_qp = npz_all['data_qp']
                    all_bins_qm = npz_all['data_qm']

                    npz_jet = np.load('./npz/'+filename_jet+'{:02d}.npz'.format(snapn))
                    jet_bins_h_dr = npz_jet['data_bin']
                    jet_bins_H_adv = npz_jet['data_H_adv']
                    jet_bins_H_conv = npz_jet['data_H_conv']
                    jet_bins_vol_approx = npz_jet['data_vol_approx']
                    jet_bins_qp = npz_jet['data_qp']
                    jet_bins_qm = npz_jet['data_qm']

                    assert (all_bins_h_dr == jet_bins_h_dr).all(), 'Segment bins of jet and shell region should be the same!'
                    amb_bins_h_dr = np.copy(jet_bins_h_dr)

                    all_adv = (np.array(all_bins_H_adv)[:,0]+np.array(all_bins_H_adv)[:,1])
                    all_conv = (np.array(all_bins_H_conv)[:,0]+np.array(all_bins_H_conv)[:,1])
                    all_qp = (np.array(all_bins_qp)[:,0]-np.array(all_bins_qp)[:,1])
                    all_qm = (np.array(all_bins_qm)[:,0]-np.array(all_bins_qm)[:,1])
                    all_dV = (np.array(all_bins_vol_approx)[:,0]-np.array(all_bins_vol_approx)[:,1])
                    all_adv_approx = all_adv / all_dV
                    all_conv_approx = all_conv / all_dV
                    all_qp_approx = all_qp / all_dV
                    all_qm_approx = all_qm / all_dV
                    all_mean_H_adv += all_adv_approx
                    all_mean_H_conv += all_conv_approx
                    all_mean_qp_approx += all_qp_approx
                    all_mean_qm_approx += all_qm_approx
                    all_mean_count += 1


                    jet_adv = (np.array(jet_bins_H_adv)[:,0]+np.array(jet_bins_H_adv)[:,1]+np.array(jet_bins_H_adv)[:,2])
                    jet_conv = (np.array(jet_bins_H_conv)[:,0]+np.array(jet_bins_H_conv)[:,1]+np.array(jet_bins_H_conv)[:,2])
                    jet_qp = (np.array(jet_bins_qp)[:,0]-np.array(jet_bins_qp)[:,1]+np.array(jet_bins_qp)[:,2])
                    jet_qm = (np.array(jet_bins_qm)[:,0]-np.array(jet_bins_qm)[:,1]+np.array(jet_bins_qm)[:,2])
                    jet_dV = (np.array(jet_bins_vol_approx)[:,0]-np.array(jet_bins_vol_approx)[:,1]+np.array(jet_bins_vol_approx)[:,2])
                    jet_adv_approx = jet_adv / jet_dV
                    jet_conv_approx = jet_conv / jet_dV
                    jet_qp_approx = jet_qp / jet_dV
                    jet_qm_approx = jet_qm / jet_dV
                    jet_mean_H_adv += jet_adv_approx
                    jet_mean_H_conv += jet_conv_approx
                    jet_mean_qp_approx += jet_qp_approx
                    jet_mean_qm_approx += jet_qm_approx
                    jet_mean_count += 1         


                    amb_adv = (all_adv-jet_adv)
                    amb_conv = (all_conv-jet_conv)
                    amb_qp = (all_qp-jet_qp)
                    amb_qm = (all_qm-jet_qm)
                    amb_dV = (all_dV-jet_dV)
                    amb_adv_approx = amb_adv / amb_dV
                    amb_conv_approx = amb_conv / amb_dV
                    amb_qp_approx = amb_qp / amb_dV
                    amb_qm_approx = amb_qm / amb_dV
                    amb_mean_H_adv += amb_adv_approx
                    amb_mean_H_conv += amb_conv_approx
                    amb_mean_qp_approx += amb_qp_approx
                    amb_mean_qm_approx += amb_qm_approx
                    amb_mean_count += 1

                    for i,adv in enumerate(amb_adv_approx):
                        if adv >= 0.:
                            amb_bins_H_adv_pos.append([amb_bins_h_dr[i], adv])
                        else:
                            amb_bins_H_adv_neg.append([amb_bins_h_dr[i], -adv])
                    for i,conv in enumerate(amb_conv_approx):
                        if conv >= 0.:
                            amb_bins_H_conv_pos.append([amb_bins_h_dr[i], conv])
                        else:
                            amb_bins_H_conv_neg.append([amb_bins_h_dr[i], -conv])

                    np.savez('./npz/'+filename_amb+'{:02d}.npz'.format(snapn), data_bin=amb_bins_h_dr, data_H_adv=amb_adv_approx, data_H_conv=amb_conv_approx, data_H_qp=amb_qp_approx, data_H_qm=amb_qm_approx, data_H_adv_pos=amb_bins_H_adv_pos, data_H_adv_neg=amb_bins_H_adv_neg, data_H_conv_pos=amb_bins_H_conv_pos, data_H_conv_neg=amb_bins_H_conv_neg)


                fig = plt.figure(num=1)
                plt.clf()
                fig, ax = plt.subplots(num=1)

                ax.plot(all_bins_h_dr, np.abs(all_qm_approx), color='C0', ls='--', label=r'$Q^-$')
                ax.plot(all_bins_h_dr, np.abs(all_qp_approx), color='C1', ls='-', label=r'$Q^+_{\mathrm{lim}}$')
                ax.plot(all_bins_h_dr, np.abs(all_adv_approx), color='C2', ls='dotted', label=r'$H_{\mathrm{adv}}$')
                ax.plot(all_bins_h_dr, np.abs(all_conv_approx), color='C3', ls='dashdot', label=r'$H_{\mathrm{conv}}$')
                
                ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
                ax.set_ylabel(r'$H_{\mathrm{all}},C_{\mathrm{all}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
                ax.set_yscale('log')
                ax.set_xlim([0,binmax])
                ax.set_ylim([1e-28,1e-22])
                ax.legend()
                ax.set_title(r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
                fig.tight_layout()
                plt.savefig('./temp_im/'+filename_all+'{:02d}.png'.format(snapn), dpi=300)
                plt.close('all')


                fig = plt.figure(num=1)
                plt.clf()
                fig, ax = plt.subplots(num=1)

                ax.plot(jet_bins_h_dr, np.abs(jet_qm_approx), color='C0', ls='--', label=r'$Q^-$')
                ax.plot(jet_bins_h_dr, np.abs(jet_qp_approx), color='C1', ls='-', label=r'$Q^+_{\mathrm{lim}}$')
                ax.plot(jet_bins_h_dr, np.abs(jet_adv_approx), color='C2', ls='dotted', label=r'$H_{\mathrm{adv}}$')
                ax.plot(jet_bins_h_dr, np.abs(jet_conv_approx), color='C3', ls='dashdot', label=r'$H_{\mathrm{conv}}$')

                ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
                ax.set_ylabel(r'$H_{\mathrm{jet}},C_{\mathrm{jet}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
                ax.set_yscale('log')
                ax.set_xlim([0,binmax])
                ax.set_ylim([1e-28,1e-22])
                ax.legend()
                ax.set_title(r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
                fig.tight_layout()
                plt.savefig('./temp_im/'+filename_jet+'{:02d}.png'.format(snapn), dpi=300)
                plt.close('all')


                fig = plt.figure(num=1)
                plt.clf()
                fig, ax = plt.subplots(num=1)

                ax.plot(amb_bins_h_dr, np.abs(amb_qm_approx), color='C0', ls='--', label=r'$Q^-$')
                ax.plot(amb_bins_h_dr, np.abs(amb_qp_approx), color='C1', ls='-', label=r'$Q^+_{\mathrm{lim}}$')

                if len(amb_bins_H_adv_pos) > 0:
                    ax.scatter(list(np.array(amb_bins_H_adv_pos)[:,0]),  list(np.array(amb_bins_H_adv_pos)[:,1]), color='r', marker='+', label=r'$H^+_{\mathrm{adv}}$')
                if len(amb_bins_H_adv_neg) > 0:
                    ax.scatter(list(np.array(amb_bins_H_adv_neg)[:,0]),  list(np.array(amb_bins_H_adv_neg)[:,1]), color='b', marker='+', label=r'$H^-_{\mathrm{adv}}$')
                if len(amb_bins_H_conv_pos) > 0:
                    ax.scatter(list(np.array(amb_bins_H_conv_pos)[:,0]), list(np.array(amb_bins_H_conv_pos)[:,1]), facecolors='none', edgecolors='r', label=r'$H^+_{\mathrm{conv}}$')
                if len(amb_bins_H_conv_neg) > 0:
                    ax.scatter(list(np.array(amb_bins_H_conv_neg)[:,0]), list(np.array(amb_bins_H_conv_neg)[:,1]), facecolors='none', edgecolors='b', label=r'$H^-_{\mathrm{conv}}$')

                ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
                ax.set_ylabel(r'$H_{\mathrm{amb}},C_{\mathrm{amb}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
                ax.set_yscale('log')
                ax.set_xlim([0,binmax])
                ax.set_ylim([1e-28,1e-22])
                ax.legend()
                ax.set_title(r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
                fig.tight_layout()
                plt.savefig('./images/'+filename_amb+'{:02d}.png'.format(snapn), dpi=300)
                plt.close('all')


                fig = plt.figure(num=1)
                plt.clf()
                fig, ax = plt.subplots(num=1)

                ax.plot(amb_bins_h_dr, np.abs(amb_qm_approx), color='C0', ls='--', label=r'$Q^-$')
                ax.plot(amb_bins_h_dr, np.abs(amb_qp_approx), color='C1', ls='-', label=r'$Q^+_{\mathrm{lim}}$')
                ax.plot(amb_bins_h_dr, np.abs(amb_adv_approx), color='C2', ls='dotted', label=r'$H_{\mathrm{adv}}$')
                ax.plot(amb_bins_h_dr, np.abs(amb_conv_approx), color='C3', ls='dashdot', label=r'$H_{\mathrm{conv}}$')

                ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
                ax.set_ylabel(r'$H_{\mathrm{amb}},C_{\mathrm{amb}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
                ax.set_yscale('log')
                ax.set_xlim([0,binmax])
                ax.set_ylim([1e-28,1e-22])
                ax.legend()
                ax.set_title(r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
                fig.tight_layout()
                plt.savefig('./temp_im/'+filename_amb+'{:02d}.png'.format(snapn), dpi=300)
                plt.close('all')

        np.savez('./npz/'+filename_all+'mean.npz', data_mean_qp_approx=all_mean_qp_approx, data_mean_qm_approx=all_mean_qm_approx, data_mean_H_adv=all_mean_H_adv, data_mean_H_conv=all_mean_H_conv, data_count=all_mean_count)
        np.savez('./npz/'+filename_jet+'mean.npz', data_mean_qp_approx=jet_mean_qp_approx, data_mean_qm_approx=jet_mean_qm_approx, data_mean_H_adv=jet_mean_H_adv, data_mean_H_conv=jet_mean_H_conv, data_count=jet_mean_count)
        np.savez('./npz/'+filename_amb+'mean.npz', data_mean_qp_approx=amb_mean_qp_approx, data_mean_qm_approx=amb_mean_qm_approx, data_mean_H_adv=amb_mean_H_adv, data_mean_H_conv=amb_mean_H_conv, data_count=amb_mean_count)


    if makemovie:
        subprocess.call('ffmpeg -y -framerate 10 -i ./temp_im/'+filename_all+'%2d.png -pix_fmt yuv420p -vf scale=1280:-2 ./temp_mv/'+filename_all+'movie.mp4', shell=True)
        subprocess.call('ffmpeg -y -i ./temp_mv/'+filename_all+'movie.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 ./temp_mv/'+filename_all+'movie.gif', shell=True)

    if makemovie:
        subprocess.call('ffmpeg -y -framerate 10 -i ./temp_im/'+filename_jet+'%2d.png -pix_fmt yuv420p -vf scale=1280:-2 ./temp_mv/'+filename_jet+'movie.mp4', shell=True)
        subprocess.call('ffmpeg -y -i ./temp_mv/'+filename_jet+'movie.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 ./temp_mv/'+filename_jet+'movie.gif', shell=True)

    if makemovie:
        subprocess.call('ffmpeg -y -framerate 10 -i ./images/'+filename_amb+'%2d.png -pix_fmt yuv420p -vf scale=1280:-2 ./movies/'+filename_amb+'movie.mp4', shell=True)
        subprocess.call('ffmpeg -y -i ./movies/'+filename_amb+'movie.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 ./movies/'+filename_amb+'movie.gif', shell=True)

        subprocess.call('ffmpeg -y -framerate 10 -i ./temp_im/'+filename_amb+'%2d.png -pix_fmt yuv420p -vf scale=1280:-2 ./temp_mv/'+filename_amb+'movie.mp4', shell=True)
        subprocess.call('ffmpeg -y -i ./temp_mv/'+filename_amb+'movie.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 ./temp_mv/'+filename_amb+'movie.gif', shell=True)


    npz_all = np.load('./npz/'+filename_all+'00.npz')
    all_bins_h_dr = npz_all['data_bin']
    npz_jet = np.load('./npz/'+filename_jet+'00.npz')
    jet_bins_h_dr = npz_jet['data_bin']
    npz_amb = np.load('./npz/'+filename_amb+'00.npz')
    amb_bins_h_dr = npz_amb['data_bin']
    assert (all_bins_h_dr == jet_bins_h_dr).all() and (jet_bins_h_dr == amb_bins_h_dr).all(), 'Segment bins of jet, shell and ambient region should be the same!'

    # just to get the time_myr:
    r_len, r_pos, time_myr, n = get_npz_data(use_h5py=False, rmax=rmax, snapnum=all_mean_count, path=path, snapvalue='n')  

    fig = plt.figure(num=1)
    plt.clf()
    fig, ax = plt.subplots(num=1)
    ax.plot(all_bins_h_dr, np.abs(all_mean_qm_approx)/all_mean_count, color='C0', ls='--', label=r'$<Q^->$')
    ax.plot(all_bins_h_dr, np.abs(all_mean_qp_approx)/all_mean_count, color='C1', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(all_bins_h_dr, np.abs(all_mean_H_adv)/all_mean_count, color='C2', ls='dotted', label=r'$<H_{\mathrm{adv}}>$')
    ax.plot(all_bins_h_dr, np.abs(all_mean_H_conv)/all_mean_count, color='C3', ls='dashdot', label=r'$<H_{\mathrm{conv}}>$')
    ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
    ax.set_ylabel(r'$H_{\mathrm{all}},C_{\mathrm{all}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
    ax.set_yscale('log')
    ax.set_xlim([0,binmax])
    ax.set_ylim([1e-28,1e-22])
    ax.legend()
    ax.set_title('mean over time '+r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
    fig.tight_layout()
    plt.savefig('./temp_im/'+filename_all+'mean.png', dpi=300)
    plt.close('all')

    fig = plt.figure(num=1)
    plt.clf()
    fig, ax = plt.subplots(num=1)
    ax.plot(jet_bins_h_dr, np.abs(jet_mean_qm_approx)/jet_mean_count, color='C0', ls='--', label=r'$<Q^->$')
    ax.plot(jet_bins_h_dr, np.abs(jet_mean_qp_approx)/jet_mean_count, color='C1', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(jet_bins_h_dr, np.abs(jet_mean_H_adv)/jet_mean_count, color='C2', ls='dotted', label=r'$<H_{\mathrm{adv}}>$')
    ax.plot(jet_bins_h_dr, np.abs(jet_mean_H_conv)/jet_mean_count, color='C3', ls='dashdot', label=r'$<H_{\mathrm{conv}}>$')
    ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
    ax.set_ylabel(r'$H_{\mathrm{jet}},C_{\mathrm{jet}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
    ax.set_yscale('log')
    ax.set_xlim([0,binmax])
    ax.set_ylim([1e-28,1e-22])
    ax.legend()
    ax.set_title('mean over time '+r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
    fig.tight_layout()
    plt.savefig('./temp_im/'+filename_jet+'mean.png', dpi=300)
    plt.close('all')

    fig = plt.figure(num=1)
    plt.clf()
    fig, ax = plt.subplots(num=1)
    ax.plot(amb_bins_h_dr, np.abs(amb_mean_qm_approx)/amb_mean_count, color='C0', ls='--', label=r'$<Q^->$')
    ax.plot(amb_bins_h_dr, np.abs(amb_mean_qp_approx)/amb_mean_count, color='C1', ls='-', label=r'$<Q^+_{\mathrm{lim}}>$')
    ax.plot(amb_bins_h_dr, np.abs(amb_mean_H_adv)/amb_mean_count, color='C2', ls='dotted', label=r'$<H_{\mathrm{adv}}>$')
    ax.plot(amb_bins_h_dr, np.abs(amb_mean_H_conv)/amb_mean_count, color='C3', ls='dashdot', label=r'$<H_{\mathrm{conv}}>$')
    ax.set_xlabel(r'$r\,[\mathrm{kpc}]$')
    ax.set_ylabel(r'$H_{\mathrm{amb}},C_{\mathrm{amb}}\ [\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-3}]$')
    ax.set_yscale('log')
    ax.set_xlim([0,binmax])
    ax.set_ylim([1e-28,1e-22])
    ax.legend()
    ax.set_title('mean over time '+r'$t={:02.2f}\,\mathrm{{Myr}}$'.format(time_myr))
    fig.tight_layout()
    plt.savefig('./images/'+filename_amb+'mean.png', dpi=300)
    plt.close('all')


if __name__ == "__main__":

# Surface Region with interpolated velocity fluxes:
    #calc_surf_int_flux_region()

# actual fluxes as in Yang+Reynolds2016:
    ###calc_Q_adv_conv_sphere()
    #calc_L_adv_conv_jet_region()
    #calc_L_adv_conv_shell_region()
    calc_L_adv_conv_ambient_region()
