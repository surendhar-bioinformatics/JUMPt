import torch as tc, numpy as np, control, numpy.linalg as nl

def gamma_setter(obj,idx,x):
    obj.gamma[idx] = tc.DoubleTensor(x).to(obj.gamma.device)**.5

def createObjFunction( etaP, LysConc, ThetaF, t, ThetaT, mapper=lambda x: x ):
    """
    Generate a callable objective function that maps gamma to error.

    t is the array of timepoints
    ThetaT is the array of ThetaP and ThetaL organized as
       [ThetaP1(0) ThetaP2(0) ... ThetaPn(0) ThetaL(0)]
       [ThetaP1(t1) ThetaP2(t1) ... ThetaPn(t1) ThetaL(t1)]
       ...
       [ThetaP1(tm) ThetaP2(tm) ... ThetaPn(tm) ThetaL(tm)]
    """
    n = ThetaT.shape[1] - 1
    ThetaT = tc.transpose(tc.DoubleTensor(ThetaT),0,1)
    ThetaTGPU = ThetaT
    t = tc.DoubleTensor(t)
    return lambda c: ((mapper(coupledOde( n, etaP, LysConc, ThetaT[:-1,0], np.array([ThetaT[-1,0]]),
                                   np.array([ThetaF]), c ).forward(t)) - mapper(ThetaTGPU))**6).sum().item()*1e3

class coupledOde:
    def __init__( self, nProteins, etaP, LysConc, ThetaP0, ThetaL0, ThetaF0, initGamma ):
        super(coupledOde,self).__init__()
        #self.gamma = (initGamma)**.5
        self.gamma = abs(initGamma)**.5
        #self.gamma = tc.nn.Parameter(tc.abs(tc.DoubleTensor(initGamma))**.5)
        self.gammaPidx = np.arange(nProteins)
        self.gammaLidx = nProteins
        self.c = etaP/LysConc
        self.theta0 = np.array(np.vstack((ThetaP0.reshape((-1,1)),ThetaL0.reshape((-1,1)),
                                          ThetaF0.reshape((-1,1)))))
        
    def A( self ):
        A = np.zeros((self.gamma.shape[0]+1,self.gamma.shape[0]+1))
        A[range(A.shape[0]-2),range(A.shape[1]-2)] = -(self.gamma[:-1]**2)
        A[:-2,-2] = self.gamma[:-1]**2
        A[-2,:-2] = self.c*(self.gamma[:-1]**2)
        A[-2,-2] = -(self.gamma[-1]**2) - (self.c*(self.gamma[:-1]**2)).sum()
        A[-2,-1] = self.gamma[-1]**2
        return A
    
    def integrate( self, t ):
        Lambda,U = nl.eig(-self.A())
        
        x = np.dot(nl.inv(U),self.theta0)
        y = x*np.exp(-t*Lambda.reshape((-1,1)))
        Z = np.dot(U,y)

        return Z[:-1,:]
    
def generateUinv( Ureal, Uimag ):
    UinvReal = tc.pinverse(Ureal + tc.mm(tc.mm(Uimag,tc.pinverse(Ureal)),Uimag),rcond=1e-9)
    UinvImag = tc.pinverse(Uimag + tc.mm(tc.mm(Ureal,tc.pinverse(Uimag)),Ureal),rcond=1e-9)    
    return UinvReal,UinvImag
    
def generateU( Lambda, U ):
    i = 0
    Ureal = tc.DoubleTensor().new_zeros(U.shape,device=U.device)
    Uimag = tc.DoubleTensor().new_zeros(U.shape,device=U.device)
    while i < Lambda.shape[0]:
        if i+1 < Lambda.shape[0] and tc.abs(Lambda[i,0] - Lambda[i+1,0]) < 1e-17 and tc.abs(Lambda[i,1] + Lambda[i+1,1]) < 1e-17:
            Ureal[:,i] = Ureal[:,i+1] = U[:,i]
            Uimag[:,i] = U[:,i+1]
            Uimag[:,i+1] = -U[:,i+1]
            i += 2
        else:
            Ureal[:,i] = U[:,i]
            i += 1

    return (Ureal,Uimag)

def cmul( areal, aimag, breal, bimag, mop ):
    return (mop(areal,breal),
            mop(areal,bimag) + mop(aimag,breal))

def cdiv( areal, aimag, breal, bimag, dop, mop ):
    return (dop(areal,breal.reshape((1,-1))),mop(aimag,0))
