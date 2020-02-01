import torch as tc, numpy as np

class MultiExp(tc.nn.Module):
    def __init__( self, A, b, xinit ):
        super(MultiExp,self).__init__()
        self.A = A
        self.b = b
        x = tc.DoubleTensor(xinit)
        x = tc.nn.Parameter(x)
        self.register_parameter('x',x)
        
    def forward( self, t ):
        return tc.mm(self.A,(self.x*(self.x > 0).to(tc.double)) + (self.x < 0).to(tc.double)*1e-4)

def fs( Lambda, C, s, backend=np ):
    return ((Lambda.reshape(list(Lambda.shape)+[1])+s.reshape((1,1,-1)))**-1*C.reshape(list(C.shape)+[1])).sum(0)
    
def ft( Lambda, C, t ):
    return np.multiply(np.exp(-Lambda.reshape(list(Lambda.shape)+[1])*t.reshape((1,1,-1))),C.reshape(list(C.shape)+[1])).sum(0)

def dft( Lambda, C, t ):
    return np.multiply(np.multiply(-Lambda.reshape(list(Lambda.shape)+[1]),
                                   np.exp(-Lambda.reshape(list(Lambda.shape)+[1])*t.reshape((1,1,-1)))),
                       C.reshape(list(C.shape)+[1])).sum(0)
