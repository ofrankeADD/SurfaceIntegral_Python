import numpy as np

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