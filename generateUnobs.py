import multiExp as me, numpy as np, control as co

class objectiveFunction:
    def __init__( self, tP, fP, tfL, dfL ):
        self.t = tP
        self.fP = fP
        self.dfL = dfL
        self.tfL = tfL
        self.mask = np.array([(tfL == tP[i]).nonzero()[0] for i in range(tP.shape[0])])

    def __call__( self, gamma ):
        residual = self.fP.reshape((-1,)) - generateSys( self.tfL, self.dfL, gamma )[self.mask].reshape((-1,))
        return 1e3*(residual**6).sum()

def generateSys( t, dfL, gamma ):
    A = np.zeros((2,2))
    A[0,0] = -gamma
    A[0,1] = gamma
    B = np.zeros((2,2))
    B[1,1] = 1
    C = np.eye(2)
    D = np.zeros((2,2))
    csys = co.ss( A, B, C, D )
    x0 = np.ones((2,1))
    U = np.vstack((np.zeros((1,dfL.shape[0])),dfL))
    rv = co.forced_response( csys, t, X0=x0, U=U )
    return rv[1][0,:]

