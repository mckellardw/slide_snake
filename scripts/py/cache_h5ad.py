import sys
from scanpy import read_10x_mtx

args = sys.argv
mat_in = args[1]
ad_out = args[2]
var_names = args[3]

adata = read_10x_mtx(
    path=mat_in,
    var_names=var_names,
    make_unique=True,
    cache=False
)

adata.write(ad_out)
