# pack-mm.py

from scipy.io import mmread

adj = mmread('data/jhu.mtx').tocsr()

n_nodes = adj.shape[0]
n_edges = adj.nnz

packed_data = np.hstack([
    np.int64(n_nodes),
    np.int64(n_edges),
    adj.indptr.astype(np.int64),
    adj.indices.astype(np.int64),
    adj.data.astype(np.int64)
])

_  = open('jhu.bin', 'wb').write(bytearray(packed_data))
