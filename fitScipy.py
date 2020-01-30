import numpy.linalg as nl, numpy as np, numpy.random as nr, parseH as ph, scipy.optimize as so, sys, h5py, tempfile, os

class NoFun:
    def __init__( self, t, n ):
        self.t = t
        self.cidx = np.arange(0,n,2)
        self.lidx = np.arange(1,n,2)
        self.lmask = np.arange(0,n) % 2
        self.cmask = self.lmask == 0
        self.expander = (np.arange(n)//2)*2
    
    def __call__( *args ):
        self = args[0]
        z = args[1]
        if len(args) > 2:
            t = args[2].reshape((1,-1))
        else:
            t = self.t
        return (np.abs(z[self.cidx]).reshape((-1,1))*np.exp(-np.abs(z[self.lidx]).reshape((-1,1))*t)).sum(0)
    
    def grad( self, x, w, t=None ):
        if None == t:
            t = self.t
        x = x.reshape((-1,1))
        return ((2*w**2*self(x[:,0]))*(np.exp(x[self.expander+1]*t.reshape((1,-1)))*(x[self.expander].T*self.lmask).T*(t*self.lmask.reshape((-1,1))))).sum(0)
            
def get_res( nofun, morefun, nexp, xhat, ft, constraintT, sign=1., rank=0 ):
    w = np.ones(ft.shape)
    w[0] *= 10.
    print( 'rank {}'.format(rank), end='' )
    res = so.minimize( lambda z: nl.norm(np.multiply(w,nofun(z) - ft))**2, xhat, bounds=[(None,None),(2**-11,None)]*nexp, jac='2-point', hess=so.BFGS(), method='trust-constr', options={'verbose':1,'maxiter':2048,'xtol':1e-8} )
    print( 'rank {}'.format(rank), end='' )    
    res = so.minimize( lambda z: nl.norm(np.multiply(w,nofun(z) - ft))**2, xhat, constraints=so.NonlinearConstraint(lambda z: sign*(nofun(z,constraintT) - morefun(constraintT)),0,np.inf), bounds=[(None,None),(2**-11,None)]*nexp, jac='2-point', hess=so.BFGS(), method='trust-constr', options={'verbose':1,'maxiter':2048,'xtol':1e-8} )
    i = 1
    while res.status == 3 or res.fun > 5e-1:
        print( 'rank {}'.format(rank), end='' )
        alpha = 1./np.sqrt(++i)
        res = so.minimize( lambda z: nl.norm(np.multiply(w,nofun(z) - ft))**2, res.x*(1. - alpha) + alpha*nr.uniform(-1e-1,1e-1,xhat.shape[0]), constraints=so.NonlinearConstraint(lambda z: sign*(nofun(z,constraintT) - morefun(constraintT)),0,np.inf), bounds=[(None,None),(2**-11,None)]*nexp, jac='2-point', hess=so.BFGS(), method='trust-constr', options={'verbose':1,'maxiter':2048,'xtol':1e-8} )
    return res

heterogeneous = False
while '--heterogeneous' in sys.argv:
    heterogeneous = True
    sys.argv.pop(sys.argv.index('--heterogeneous'))

nexp = 3
only = None
for i in range(len(sys.argv)):
    if '--nexp' in sys.argv[i]:
        k = sys.argv.pop(i)
        nexp = int(k.split('=')[1])
        break
for i in range(len(sys.argv)):    
    if '--only' in sys.argv[i]:
        k = sys.argv.pop(i)
        only = k.split('=')[1]
        break
        
    
if heterogeneous == True:
    t,ft,header = ph.parseH(sys.argv[1])
else: 
    reader = csv.reader(open(sys.argv[1]))
    lines = list(reader)
    header = lines[0]
    data = np.array(lines[1:]).astype(np.double)
    ft = [data[:,i] for i in range(1,data.shape[1])]
    t = [data[:,0] for i in range(len(ft))]


resl = []
w = np.ones(ft[0].shape)
w[0] *= 10.
mins = dict(zip(t[0],ft[0]))
for time,ftime in zip(t[1:],ft[1:]):
    for i in range(len(time)):
        mins[time[i]] = min(ftime[i],mins[time[i]])

nr.seed(0)
nofun = lambda z,t=t[0].reshape((1,-1)): (np.abs(z[np.arange(0,z.shape[0],2)].reshape((-1,1)))*np.exp(-np.abs(z[np.arange(1,z.shape[0],2)]).reshape((-1,1))*t)).sum(0)
morefun = lambda x: np.array([mins[t] for t in sorted(mins.keys())])
xhat = nr.uniform(0,.5,nexp*2)
if None != only:
    res = get_res( nofun, morefun, nexp, xhat, ft[0], np.array(list(mins.keys())), sign=-1., rank=only )
    for j in range(8):
        xhat = nr.uniform(0,.5,nexp*2)
        trial = get_res( nofun, morefun, nexp, xhat, ft[0], np.array(list(mins.keys())), sign=-1., rank=only )
        if trial.fun < res.fun:
            res = trial
        if res.fun < 1e-5:
            break
    resl.append( res )
    td = nr.uniform(min(mins.keys()),max(mins.keys()),128)
    i = header.index(only)
    nofun = NoFun( t[i], nexp*2 )
    morefun = lambda t: (np.abs(resl[0].x[np.arange(0,resl[0].x.shape[0],2)]).reshape((-1,1))*np.exp(-np.abs(resl[0].x[np.arange(1,resl[0].x.shape[0],2)]).reshape((-1,1))*t)).sum(0)
    xhat = nr.uniform(0,.5,resl[0].x.shape[0])
    res = get_res( nofun, morefun, nexp, xhat, ft[i], td, rank=only )

    for j in range(8):
        xhat = res.x + nr.uniform(-1e-3,1e-3,resl[0].x.shape[0])
        trial = get_res( nofun, morefun, nexp, xhat, ft[i], td, rank=only )
        if trial.fun < res.fun:
            res = trial
        if res.fun < 1e-5:
            break
        
    resl.append( res )
    of = h5py.File( sys.argv[2], 'w' )
    of.create_dataset( 'C', data=np.abs(np.hstack([res.x[np.arange(0,res.x.shape[0],2)].reshape((-1,1)) for res in resl])) )
    of.create_dataset( 'Lambda', data=np.abs(np.hstack([res.x[np.arange(1,res.x.shape[0],2)].reshape((-1,1)) for res in resl])) )
    of.create_dataset( 'td', data=td )    
    of.close()
else:
    mf = tempfile.NamedTemporaryFile(mode='w')
    mf.file.write( 'all: {}\n'.format(str.join(' ',[pid+'.h5' for pid in header[1:]])) )
    for pid in header[1:]:
        mf.file.write( '{}:\n\tpython fitScipy.py --nexp={} {} --only={} {} {}\n'.format(pid+'.h5',nexp,
                                                                            '--heterogeneous' if heterogeneous
                                                                            else '',
                                                                            pid,sys.argv[1],pid+'.h5') )
    mf.file.flush()
    os.system( 'make -f {} -j 20'.format(mf.name) )
    hf = h5py.File( header[1]+'.h5', 'r' )
    Lambda = [np.array(hf['Lambda'])[:,:2]]
    C = [np.array(hf['C'])[:,:2]]
    for pid in header[2:]:
        hf = h5py.File( pid+'.h5', 'r' )
        Lambda.append(np.array(hf['Lambda'])[:,1].reshape((-1,1)))
        C.append(np.array(hf['C'])[:,1].reshape((-1,1)))

    of = h5py.File( sys.argv[2], 'w' )
    of.create_dataset( 'td', data=hf['td'] )
    of.create_dataset( 'C', data=np.hstack(C) )
    of.create_dataset( 'Lambda', data=np.hstack(Lambda) )
    of.close()
        
