import csv, re, numpy as np

def squeeze( col ):
    idx = (1 - np.isnan(col)).nonzero()
    return col[idx]

def parseH( fname, do_squeeze=True ):
    reader = csv.reader(open(fname))
    l = [l for l in reader]
    header = l[0]
    fx = np.array([t[1:] for t in l if re.search('^Day',t[0])]).astype(np.double)
    x = np.array([t[1:] for t in l if re.search('^T_Day',t[0])]).astype(np.double)
    if do_squeeze:
        t = [squeeze(x[:,i]) for i in range(x.shape[1])]
        ft = [squeeze(fx[:,i]) for i in range(fx.shape[1])]
    else:
        t = x
        ft = fx
        
    return t,ft,header[1:]
