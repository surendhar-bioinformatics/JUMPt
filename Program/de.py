import numpy as np, scipy.sparse as ss, numpy.linalg as nl, scipy.fftpack as sft, torch as tc

def linSys( ThetaHat, ThetaA, dThetaHat, dThetaA, etaP_LysC ):
    """
    Generate one system for a time point.

    ThetaHat := (Theta_L(t) - Theta_p(t))  
    ThetaHat is num_protiens x 1
    ThetaA is a scalar := (Theta_F(t) - Theta_L(t))
    etaP_LysC is num_protiens x 1, contains the eta_P/Lysine concentration ratio

    output: A,b for min_x \|Ax - b\|_2^2
    """
    n = ThetaHat.shape[0] + 1
    r = np.hstack((np.arange(n),
                   (n-1)*np.ones(n-1)))
    c = np.hstack((np.arange(n),
                   np.arange(n-1)))
    v = np.hstack((ThetaHat.T,np.ones(1)*ThetaA,-ThetaHat.T*etaP_LysC))
    A = ss.coo_matrix( (v,(r,c)), shape=(n,n) )
    b = np.hstack((dThetaHat.T,np.ones([1]*len(dThetaHat.shape))*dThetaA))
    return (A,np.matrix(b).T)

def lppBoltOn( n, blocksize, t, factor=.05 ):
    r = (np.arange(1,n+1)*blocksize) - 1
    c = np.arange(n)
    v = np.exp(t*factor)
    return ss.coo_matrix( (v,(r,c)), shape=(n*blocksize,n) ).tocsr()

def lstSqMultiPt( ThetaHat, ThetaA, dThetaHat, dThetaA, etaP_LysC, t, factor=.05 ):
    """
    Generate one system for m time points.

    ThetaHat := (Theta_L(t) - Theta_p(t))  
    ThetaHat is num_protiens x m (where m is the number of time points)
    ThetaA is m x 1  := (Theta_F(t) - Theta_L(t))
    etaP_LysC is num_protiens x 1, contains the eta_P/Lysine concentration ratio
    t is time

    output: A,b,funpack for min_x \|Ax - b\|_2^2
            funpack(x) is a function that returns (gamma_p,gamma_a) when given x.
    """
    Abl = [linSys( ThetaHat[i,:], ThetaA[i], dThetaHat[i,:], dThetaA[i], etaP_LysC )
           for i in range(ThetaHat.shape[0])]
    return (ss.vstack([tu[0]*np.exp(ti*factor) for tu,ti in zip(Abl,t)]),np.vstack([tu[1]*np.exp(ti*factor) for tu,ti in zip(Abl,t)]),
            lambda x: (x[:-1],x[-1]))

def lstSqMultiPtwPU( ThetaHat, ThetaA, dThetaHat, dThetaA, etaP_LysC, t, etaPU_LysC, lambdaPU ):
    """
    Generate one system for m time points with unknown lumped protien pool.  
    Unknown lumped protein is modeled with concetration etaPU and decay lambdaPU.  
    t is the time points on which the functions are sampled.

    ThetaHat := (Theta_L(t) - Theta_p(t))  
    ThetaHat is num_protiens x m (where m is the number of time points)
    ThetaA is m x 1  := (Theta_F(t) - Theta_L(t))
    etaP_LysC is num_protiens x 1, contains the eta_P/Lysine concentration ratio

    output: A,b,funpack for min_x \|Ax - b\|_2^2
            funpack(x) is a function that returns (gamma_p,gamma_a) when given x.
    """
    ThetaHat2 = np.hstack((ThetaHat,np.exp(lambdaPU*t.reshape((-1,1)))))
    dThetaHat2 = np.hstack((dThetaHat,lambdaPU*np.exp(lambdaPU*t.reshape((-1,1)))))
    etaP_LysC2 = np.hstack((etaP_LysC,[etaPU_LysC]))
    Abl = [linSys( ThetaHat2[i,:], ThetaA[i], dThetaHat2[i,:], dThetaA[i], etaP_LysC2 )
           for i in range(ThetaHat.shape[0])]
    return (ss.vstack([t[0] for t in Abl]).tocsc(),np.vstack([t[1] for t in Abl]),
            lambda x: (x[:-1],x[-1]))

def solveSparseLstSq( A, b, W=None ):
    """
    Solve a sparse least squares problem \min_x \|Ax - b\|_2^2

    Uses the normal equations A.TAx = b
    """
    if None == W:
        W = ss.identity(A.shape[0])        
    G = (A.T*(W*A)).todense()
    return G**-1*(A.T*(W*b))

class DiffFitter(tc.nn.Module):
    def __init__( self, nGamma, nUt ):
        super(DiffFitter,self).__init__()
        self.gamma = tc.nn.Parameter(tc.DoubleTensor().new_empty((nGamma,1)).uniform_())
        self.ut = tc.nn.Parameter(tc.DoubleTensor().new_empty((nUt,1)).uniform_())
        self.register_parameter('gamma',self.gamma)
        self.register_parameter('ut',self.ut)

    def forward( self, A, B ):
        return tc.mm(A,self.gamma) + tc.mm(B,self.ut)

        
