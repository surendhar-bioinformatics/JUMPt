import numpy as np, numpy.fft as nft, scipy.interpolate as si

def diffFourier( t, ThetaT, samp=1e6 ):
    xd = (np.arange(int(samp)*np.ceil(t[-1])) + 1)/samp
    ThetaTr = np.empty((xd.shape[0],ThetaT.shape[1]))
    
    for i in range(ThetaTr.shape[1]):
        finterp = si.interp1d( t, ThetaT[:,i], kind=2 )
        ThetaTr[:,i] = finterp( xd )

    F = np.hstack([nft.rfft( np.hstack((ThetaTr[:,i],ThetaTr[-np.arange(1,ThetaTr.shape[0]+1),i])) ).reshape((-1,1))
                   for i in range(ThetaTr.shape[1])])
    omega = (nft.rfftfreq(ThetaTr.shape[0]*2,d=samp**-1)*1j*np.pi*2).reshape((-1,1))
    F *= omega

    return (xd,ThetaTr,np.hstack([nft.irfft(F[:,i]).reshape((-1,1)) for i in range(F.shape[1])]))
     
