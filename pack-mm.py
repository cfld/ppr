# pack-mm.py

from scipy.io import mmread

inpath  = 'data/jhu.mtx'
outpath = 'data/jhu.bin'
adj     = mmread(inpath).tocsr()

n_nodes = adj.shape[0]
n_edges = adj.nnz

packed_data = np.hstack([
    np.int64(n_nodes),
    np.int64(n_edges),
    adj.indptr.astype(np.int64),
    adj.indices.astype(np.int64)
])

_  = open(outpath, 'wb').write(bytearray(packed_data))
