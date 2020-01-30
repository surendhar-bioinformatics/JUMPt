import torch as tc

def fitReg( reg, t, ft, maxiter=100, tol=1e-7, w=None ):
    optimizer = tc.optim.Adamax(reg.parameters(),lr=.001)
    criterion = tc.nn.MSELoss()

    def closure():
        optimizer.zero_grad()
        output = reg.forward(t)
        loss = criterion( output, ft )
        loss.backward()
        return loss

    U = tc.DoubleTensor().new_empty((ft.shape[0],ft.shape[0]),device=ft.device)
    for n in range(maxiter):
        optimizer.zero_grad()
        output = reg.forward(t)
        if isinstance(w,type(None)):
            w = tc.DoubleTensor().new_empty((128,1,ft.shape[1]),dtype=ft.dtype,device=ft.device).uniform_()
        loss = criterion( w*output, w*ft )
        loss.backward()
        print(reg.gamma.grad)
        optimizer.step(closure)
        print( 'n = {}, loss = {}'.format(n,loss.item()) )

        if tc.abs(tc.abs(loss)).item() < tol:
            break

    return reg


def fitDiff( reg, A, B, c, maxiter=100, tol=1e-7, w=None ):
    optimizer = tc.optim.LBFGS(reg.parameters(),lr=.1)
    criterion = tc.nn.MSELoss()

    def closure():
        optimizer.zero_grad()
        output = reg.forward(A,B)
        loss = criterion( output, c )
        loss.backward()
        return loss

    for n in range(maxiter):
        optimizer.zero_grad()
        output = reg.forward(A,B)
        loss = criterion( output, c ) + tc.nn.functional.relu(-reg.gamma).sum()
        loss.backward()
        optimizer.step(closure)
        print( 'n = {}, loss = {}'.format(n,loss.item()) )

        if tc.abs(tc.abs(loss)).item() < tol:
            break

    return reg


