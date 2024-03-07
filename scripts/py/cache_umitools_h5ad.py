import sys
from scanpy import read_umi_tools

args = sys.argv
filename = args[1]
ad_out = args[2]

adata = read_umi_tools(
    filename=filename
    # make_unique=True,
    # cache=False
)

adata.write(ad_out)
