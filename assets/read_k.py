import h5py
import pickle
import sys

hdf5file = str(sys.argv[1])

with h5py.File(hdf5file,'r') as f:

    dataset = f['datablock_0/wavepacket/basisshapes/basis_shape_-2083794085692656251/']

    int_in_bytes = dataset.attrs['K'].tobytes()
    K = pickle.loads(int_in_bytes)
    print("K =", K)

    
