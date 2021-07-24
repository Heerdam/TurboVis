import h5py
import pickle
import sys

hdf5file = str(sys.argv[1])

with h5py.File(hdf5file,'r') as f:
    print("P:")
    P_list = f['datablock_0/wavepacket/Pi/P']
    for P in P_list:
        print(P)
    print("\n")

    print("Q:")
    Q_list = f['datablock_0/wavepacket/Pi/Q']
    for Q in Q_list:
        print(Q)
    print("\n")

    print("p:")
    p_list = f['datablock_0/wavepacket/Pi/p']
    for p in p_list:
        print(p)
    print("\n")

    print("q:")
    q_list = f['datablock_0/wavepacket/Pi/q']
    for q in q_list:
        print(q)
    print("\n")

    dataset = f['datablock_0/wavepacket/basisshapes/basis_shape_-2083794085692656251/']

    int_in_bytes = dataset.attrs['K'].tobytes()
    K = pickle.loads(int_in_bytes)
    print("K =", K)

    int_in_bytes = dataset.attrs['dimension'].tobytes()
    dimension = pickle.loads(int_in_bytes)
    print("dimension =", dimension)

    string_in_bytes = dataset.attrs['type'].tobytes()
    basis_type = pickle.loads(string_in_bytes)
    print("basis_type =", basis_type)
    print("\n")

    c_list = f['datablock_0/wavepacket/coefficients/c_0']
    i = 0
    for c in c_list:
        print("Step {0}:\n".format(i), c)
        i += 1
        print("\n")

    
