import spreg
import scipy.sparse
from libpysal.weights import WSP
from libpysal.weights import WSP2W


def run_se_GMM_py(y, x, w):
    sparse_w = scipy.sparse.csr_matrix(w)
    wsp = WSP(sparse_w)
    w_reg = WSP2W(wsp)
    model = spreg.GM_Error_Het(y, x, w=w_reg)

    return model
