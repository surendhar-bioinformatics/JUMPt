import numpy as np, scipy.optimize as so, scipy.integrate as si, pandas as pd, de, numpy.linalg as nl, coupledOde as co, scipy.sparse as ss, generateUnobs as gu, matplotlib.pyplot as plt, time, subprocess, os
from datetime import datetime 
from pathlib import Path
import tkinter


def runMatlab( scriptFile ):
    """
    Run matlab on scriptFile and return when finished.
    """
    retval = subprocess.run( ['matlab', '-wait', '-nojvm','-nodisplay', '-nosplash',  '-minimize' ,'-r', scriptFile ] )
    retval.check_returncode()

def justify(a, invalid_val=0, axis=1, side='left'):    
    """
    Justifies a 2D array
    
    Parameters
    ----------
    A : ndarray
        Input array to be justified
    axis : int
        Axis along which justification is to be made
    side : str
        Direction of justification. It could be 'left', 'right', 'up', 'down'
        It should be 'left' or 'right' for axis=1 and 'up' or 'down' for axis=0.
    """

    if invalid_val is np.nan:
        mask = ~np.isnan(a)
    else:
        mask = a!=invalid_val
    justified_mask = np.sort(mask,axis=axis)
    if (side=='up') | (side=='left'):
        justified_mask = np.flip(justified_mask,axis=axis)
    out = np.full(a.shape, invalid_val) 
    if axis==1:
        out[justified_mask] = a[mask]
    else:
        out.T[justified_mask.T] = a.T[mask.T]
    return out

# ====================== Function to get the sum opf squared error======================
LSE=lambda data,model:((data-model)**2).sum()

# ====================== differential-eq-system ======================
def PT_ODE(t, y0, param, A1, EtaP):

    y1 = y0[0]; y2 = np.transpose(y0[1:]); #print('y1 = ', y1); #print('y2 = ', y2)# row vector
    gama_A1 = np.abs(param[0]); gama_P = np.abs(param[1:]); #print('gama_P = ', gama_P); print('gama_A1 = ', gama_A1)
 
    # define the odes
    dy1 = gama_A1 * (0.05 - y1) + (gama_P *(EtaP/A1)*(y2 - y1)).sum(); #print('dy1 = ', dy1) #A1= Lysin a.a.
    dy2 = gama_P *(y1 - y2); #print('dy2 = ', dy2)#proteins
    
    return np.concatenate((dy1, np.transpose(dy2)), axis=None)   
    
def PT_ODE2(y0, t, param, A1, EtaP):
    y1 = y0[0]; y2 = np.transpose(y0[1:]); #print('y2 = ', y2) #print('y3 = ', y3)# row vector
    gama_A1 = np.abs(param[0]); gama_P = np.abs(param[1:]); #print('gama_A1 = ', gama_A1) #print('gama_P = ', gama_P)
    
    # define the odes
    dy1 = gama_A1 * (0.05 - y1) + (gama_P *(EtaP/A1)*(y2 - y1)).sum(); #print('dy0 = ', dy0) #A1= Lysin a.a.
    dy2 = gama_P *(y1 - y2); #print('dy1 = ', dy1)#proteins
    
    #print('np.array([dy0, np.transpose(dy1)]) = ', np.concatenate((dy0, np.transpose(dy1)), axis=None))
    return np.concatenate((dy1, np.transpose(dy2)), axis=None)   

#3.Score Fit of System
#=========================================================
def LSE_Lys_Pro(param_0, y0, EtaP, tspan, P_ratio, A1_ratio, A1):
    
    A1_P_ratio_simu = si.solve_ivp(lambda t, y: PT_ODE(t, y, param_0, A1, EtaP), [0, 32], y0, method="RK45", t_eval=tspan,); 
    LSE_Lys = ((np.transpose(A1_ratio) - A1_P_ratio_simu.y[0,:])**2).sum() ;  
    LSE_Lys_Pro = LSE_Lys + np.nansum((np.transpose(P_ratio) - A1_P_ratio_simu.y[1:,:])**2);  
    return LSE_Lys_Pro

def LSE_Lys_Pro_withUIP(param_0, y0, EtaP_UIP, tspan, P_ratio, A1_ratio, A1):
    
    start = time.time()
    A1_P_ratio_simu = si.solve_ivp(lambda t, y: PT_ODE(t, y, param_0, A1, EtaP_UIP), [0, 32], y0, method="RK45", t_eval=tspan,); 
    end = time.time()
    
    LSE_Lys = ((np.transpose(A1_ratio) - A1_P_ratio_simu.y[0,:])**2).sum() ;  
    residual = (np.transpose(P_ratio) - A1_P_ratio_simu.y[1:-1,:]) ; 
    mask = 1 - np.isnan(residual); 
    LSE_Lys_Pro = LSE_Lys + (residual[mask]**2).sum();  
    end2 = time.time()
    print('Time for ODE (incl UIP) excution (RK45) is = {:.3f} ' .format(end - start), 'and residual func exc. is = {:.3f}' .format(end2 - end), 'sec and resi =', LSE_Lys_Pro)
    return LSE_Lys_Pro

   
################################################################### 
#================== 1.Load the original data ======================
################################################################### 
# Ask the user to select a single file name.
Inputfile = input('Enter your input data file  (if it doesnt exist in current folder, enter with complete path;  test data file name is "test_data.xlsx"):') 
#root = tkinter.Tk()
#my_filetypes = [('all files', '.xlsx'), ('Excel file', '.txt')] # Build a list of tuples for each file type the file dialog should display
#Inputfile = tkinter.filedialog.askopenfilename(initialdir=os.getcwd(), title="Please select a file:", filetypes=my_filetypes)   
#root.destroy()
#Inputfile   = "test_data.xlsx" # Assign spreadsheet filename to `file`
P_data      = pd.read_excel(Inputfile,sheet_name= 'data') # Load spreadsheet
tspan       = (P_data.loc[0, ['data_1', 'data_2', 'data_3', 'data_4', 'data_5']]).values.astype(int); 
A1_ratio    = (np.transpose((P_data.loc[1, ['data_1', 'data_2', 'data_3', 'data_4', 'data_5']]).values.astype(float))).reshape((-1,1)); 
P_ratio     = np.transpose((P_data.loc[2:, ['data_1', 'data_2', 'data_3', 'data_4', 'data_5']]).values); 
EtaP        = np.transpose((P_data.loc[2:,['Lys_Conc_microM']]).values);
EtaP        = EtaP[0,:].reshape(-1)
sum_EtaP    = 183190; #micoM
A1 = 206    #microM concentration
len_P_data  = len(P_ratio[0])# len(A) is for length of rows and the number of columns is len(A[0])
tspan_medium = np.linspace(0,32,321)
tspan_long   = np.linspace(0,32,3201)      
tspan_data  = np.transpose(np.tile(tspan, (len_P_data,1 )))
ind         = ~np.isfinite(P_ratio); #ind = np.isnan(np.asarray(P_ratio.astype(float)))
tspan_data[ind] = -1

#if not os.path.exists('Results'):
#   os.makedirs('Results')
Path('Results').mkdir(parents=True, exist_ok=True)   

writer  = pd.ExcelWriter('Results/Original_data.xlsx', engine='xlsxwriter'); # Create a Pandas Excel writer using XlsxWriter as the engine.
P_data.to_excel(writer, sheet_name='data', index=False) # Convert the dataframe to an XlsxWriter Excel object.
writer.save()# Close the Pandas Excel writer and output the Excel file.
writer.close()

#=========== plot the Original data ======================
plt.figure(1)
plt.plot(tspan, A1_ratio,     'r-s',linewidth=2, markeredgecolor='r',markeredgewidth=2, markersize=10,label='A_L- (Lys unbound)')
plt.plot(justify(tspan_data[:,1], invalid_val=-1, axis=0, side='left'),justify(P_ratio[:,1],invalid_val=np.nan, axis=0, side='left'), 'b-o',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=10,label='P_L- (Lys bound protein)')
plt.plot(justify(tspan_data, invalid_val=-1, axis=0, side='left'),justify(P_ratio,invalid_val=np.nan, axis=0, side='left'), 'b-o',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=5)
plt.plot(tspan, A1_ratio, 'r-s',linewidth=2, markeredgecolor='r',markeredgewidth=2, markersize=10)
plt.title('Original data after outliers removed',FontSize = 16)
plt.xlabel('Time (days)', FontSize = 14); plt.ylabel('Fraction of  \'light\'  Lys ', FontSize=14)
plt.legend(loc='best'); plt.ylim(0, 1.03); #plt.show()
plt.savefig('Results/Original_data.jpeg', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', papertype='None',transparent='False', pad_inches=0.1, frameon='None', metadata='None', bbox_inches='tight')


################################################################## 
#==================== Individual protein optimization =============
################################################################### 
print('\n\n Calling MATLAB program to fit the indivual proteins and calculating its derivatives')
runMatlab( 'PT_indi_Fitting' )
print('\n\n Completed indiviual fitting and results were saved')
      

################################################################### 
###====== Parmeter optimization using Derivative fitting ======##
################################################################### 

Inputfile       = "Results/Indi_fit_res.xlsx" # Assign spreadsheet filename to `file`
prot_SimuTraj   = pd.read_excel(Inputfile,sheet_name= 'f_Lys_P') 
deri_prot       = pd.read_excel(Inputfile,sheet_name= 'df_Lys_P') 
f               = (prot_SimuTraj.iloc[1:,:]).values
df              = (deri_prot.iloc[1:,:]).values
t               = tspan_medium[1:] 


ThetaH = np.hstack([f[:,0].reshape((-1,1)) - f[:,i].reshape((-1,1)) 
                    for i in range(1,f.shape[1])])
ThetaA = .05-f[:,0]

A,b = de.lstSqMultiPt( ThetaH[:,:], ThetaA.reshape((-1,1)), df[:,1:], (df[:,0]-.05).reshape((-1,1)), EtaP/206, t, .0125 )[:2]
B   = de.lppBoltOn(ThetaH.shape[0],A.shape[1],t,0)
M   = ss.hstack((A,B)).tocsr()
# preconditioning with Jacobi preconditioner
D   = ss.dia_matrix( (1./np.power(M.multiply(M).sum(0),.5),0), (M.shape[1],M.shape[1]) )
G   = np.array((M.T*M).todense())
# Tikhonov matrix for regularization
MTb = np.array(M.T*b).reshape((-1,))
nidx= np.nonzero(np.dot(nl.inv(G),MTb) < 0)[0]
T   = ss.dia_matrix( (np.hstack((np.zeros(A.shape[1]),np.exp(.0125*t) -1)),0), G.shape  )
gammaUt2 = np.array(np.dot(nl.inv(G + T*1e-6), MTb)).reshape((-1,))
gammaUt = np.array(np.dot(nl.inv(G + T*1e-6), MTb)).reshape((-1,))
res     = so.minimize( lambda x: nl.norm(np.dot(G,x) - MTb.T,2)**2 + nl.norm(x*(x < 0),2)**2, np.array(np.abs(gammaUt)).reshape((-1,)), jac=lambda x: (np.dot(G,x) - MTb.T) + (x*(x < 0)), method='tnc', options={'maxiter':126,'disp':3} )
res     = so.minimize( lambda x: nl.norm(np.dot(G,x) - MTb.T,2)**2, res.x, jac=lambda x: (np.dot(G,x) - MTb.T), method='trust-constr', options={'maxiter':126,'verbose':3}, bounds=[(0,None)]*G.shape[0] )
gammaUt = res.x.reshape((-1,1))
gammaUt_temp = res.x.reshape((-1,1))
ut      = np.array(gammaUt[A.shape[1]:]).reshape((-1,1))
ut[0] = 0

def noFun( gamma ):
    fP = gu.generateSys( t, df[:,0], gamma )
    
    return np.sin(np.arccos((np.dot(gamma*(f[:,0] - fP),ut)/(nl.norm(gamma*(f[:,0] - fP))*nl.norm(ut)))))

res = so.minimize_scalar( noFun, bounds=[0,20], 
                          method='bounded', options={'print':3,'disp':3,'maxiter':2000} )
gammaUnk = res.x
ThetaUnk = gu.generateSys( t, df[:,0], gammaUnk )
dThetaUnk= gammaUnk*(f[:,0] - ThetaUnk)
c        = np.dot(ut.T,-dThetaUnk)/(nl.norm(dThetaUnk)**2)*A1
EtaP_forDeriFit = np.hstack((EtaP, c))

gamma     = np.hstack((gammaUt[:A.shape[1]].reshape((-1,)),[0]))
gamma[-1] = gamma[-2]
gamma[-2] = gammaUnk

integrator = co.coupledOde( len(EtaP_forDeriFit), EtaP_forDeriFit, A1,
                            np.ones(len(EtaP_forDeriFit)), np.ones(1), np.ones(1)*.05, gamma )


refft   = np.hstack((A1_ratio, P_ratio[:,:])); 
reft = np.hstack((np.transpose(tspan.reshape(1,-1)), tspan_data[:,:]))
reft    = np.vstack([x.reshape((1,-1)) for x in reft])
refft   = np.vstack([x.reshape((1,-1)) for x in refft])
rv      = integrator.integrate(reft[:,0])
idx     = list(range(len(EtaP_forDeriFit)-1)) + [-1]                 
residual= rv[idx,:] - refft[:,list(range(1,len(EtaP_forDeriFit)))+[0]].T
mask    = (1 - np.isnan(residual)).nonzero()
flatResidual = residual[mask]; ResError = nl.norm(flatResidual)**2;
ResError_temp = 1

while ResError_temp < ResError:
    res     = so.minimize( lambda x: nl.norm(np.dot(G,x) - MTb.T,2)**2 + nl.norm(x*(x < 0),2)**2, np.array(np.abs(gammaUt_temp)).reshape((-1,)), jac=lambda x: (np.dot(G,x) - MTb.T) + (x*(x < 0)), method='tnc', options={'maxiter':126} )
    res     = so.minimize( lambda x: nl.norm(np.dot(G,x) - MTb.T,2)**2, res.x, jac=lambda x: (np.dot(G,x) - MTb.T), method='trust-constr', options={'maxiter':126,'verbose':3}, bounds=[(0,None)]*G.shape[0] )
    gammaUt_temp = res.x.reshape((-1,1))

    ut      = np.array(gammaUt_temp[A.shape[1]:]).reshape((-1,1))
    ut[0] = 0
    
    def noFun( gamma ):
        fP = gu.generateSys( t, df[:,0], gamma )
        
        return np.sin(np.arccos((np.dot(gamma*(f[:,0] - fP),ut)/(nl.norm(gamma*(f[:,0] - fP))*nl.norm(ut)))))
    
    res = so.minimize_scalar( noFun, bounds=[0,20], 
                              method='bounded', options={'print':3,'disp':3,'maxiter':2000} )
    gammaUnk = res.x
    ThetaUnk = gu.generateSys( t, df[:,0], gammaUnk )
    dThetaUnk= gammaUnk*(f[:,0] - ThetaUnk)
    c        = np.dot(ut.T,-dThetaUnk)/(nl.norm(dThetaUnk)**2)*A1
    EtaP_forDeriFit = np.hstack((EtaP, c))
    
    gamma     = np.hstack((gammaUt_temp[:A.shape[1]].reshape((-1,)),[0]))
    gamma[-1] = gamma[-2]
    gamma[-2] = gammaUnk
    
    integrator = co.coupledOde( len(EtaP_forDeriFit), EtaP_forDeriFit, A1,
                                np.ones(len(EtaP_forDeriFit)), np.ones(1), np.ones(1)*.05, gamma )
    
    
    refft   = np.hstack((A1_ratio, P_ratio[:,:])); 
    reft = np.hstack((np.transpose(tspan.reshape(1,-1)), tspan_data[:,:]))
    reft    = np.vstack([x.reshape((1,-1)) for x in reft])
    refft   = np.vstack([x.reshape((1,-1)) for x in refft])
    rv      = integrator.integrate(reft[:,0])
    idx     = list(range(len(EtaP_forDeriFit)-1)) + [-1]                 
    residual= rv[idx,:] - refft[:,list(range(1,len(EtaP_forDeriFit)))+[0]].T
    mask    = (1 - np.isnan(residual)).nonzero()
    flatResidual = residual[mask]; ResError_temp = nl.norm(flatResidual)**2; 
    print('\n ResError_temp = ', ResError_temp,'ResError = ', ResError)
    if ResError_temp < ResError:
        ResError = ResError_temp
        gammaUt = gammaUt_temp



ut      = np.array(gammaUt[A.shape[1]:]).reshape((-1,1))
ut[0] = 0

def noFun( gamma ):
    fP = gu.generateSys( t, df[:,0], gamma )
    
    return np.sin(np.arccos((np.dot(gamma*(f[:,0] - fP),ut)/(nl.norm(gamma*(f[:,0] - fP))*nl.norm(ut)))))

res = so.minimize_scalar( noFun, bounds=[0,20], 
                          method='bounded', options={'print':3,'disp':3,'maxiter':2000} )
gammaUnk = res.x
ThetaUnk = gu.generateSys( t, df[:,0], gammaUnk )
dThetaUnk= gammaUnk*(f[:,0] - ThetaUnk)
c        = np.dot(ut.T,-dThetaUnk)/(nl.norm(dThetaUnk)**2)*A1
EtaP_forDeriFit = np.hstack((EtaP, c))

gamma     = np.hstack((gammaUt[:A.shape[1]].reshape((-1,)),[0]))
gamma[-1] = gamma[-2]
gamma[-2] = gammaUnk

integrator = co.coupledOde( len(EtaP_forDeriFit), EtaP_forDeriFit, A1,
                            np.ones(len(EtaP_forDeriFit)), np.ones(1), np.ones(1)*.05, gamma )


refft   = np.hstack((A1_ratio, P_ratio[:,:])); 
reft = np.hstack((np.transpose(tspan.reshape(1,-1)), tspan_data[:,:]))
reft    = np.vstack([x.reshape((1,-1)) for x in reft])
refft   = np.vstack([x.reshape((1,-1)) for x in refft])
rv      = integrator.integrate(reft[:,0])
idx     = list(range(len(EtaP_forDeriFit)-1)) + [-1]                 
residual= rv[idx,:] - refft[:,list(range(1,len(EtaP_forDeriFit)))+[0]].T
mask    = (1 - np.isnan(residual)).nonzero()
flatResidual = residual[mask]; ResError = nl.norm(flatResidual)**2;
print('Residual Error (LSE) after derivative fitting = ', ResError)


gamma2  = np.zeros([len(gamma)]); gamma2[0] = gamma[-1]; gamma2[1:len(gamma)] = gamma[0:-1]; gamma2 =abs(gamma2)
row     = ['gamma_P'+str(i) for i in range(1, len_P_data+1)]; row =  ['gamma_Lys'] + row +['gamma_UIP'];
df_deriFit = pd.DataFrame(gamma2,index = row, columns = ['gamma_deriFit']); 
writer  = pd.ExcelWriter('Results/Derivative_fit_res.xlsx', engine='xlsxwriter'); # Create a Pandas Excel writer using XlsxWriter as the engine.
df_deriFit.to_excel(writer, sheet_name='gamma_DeriFit') # Convert the dataframe to an XlsxWriter Excel object.
writer.save()# Close the Pandas Excel writer and output the Excel file.
writer.close()

################################################################### 
###============== Global optimization for degrdation rates======##
################################################################### 

print('\n Calling MATLAB program to fit the Lys and all proteins Globally ...')

runMatlab( 'PT_Glob_Fitting' )

print('\n Completed Global fitting and results were saved')



#########################################################################################################################  
###============ Generating figure with Global optimized parameters  and saving the tragectries of all proteins  ======##
#########################################################################################################################  

Inputfile   = Path("Results\Glob_fit_res.xlsx") # Assign spreadsheet filename to `file`
A1_P_ratio_simu  = pd.read_excel(Inputfile,sheet_name= 'traj_Lys_P') # Load spreadsheet
A1_P_ratio_simu  = A1_P_ratio_simu.values
plt.figure(1)
plt.plot(tspan_medium,A1_P_ratio_simu[:,0], 'r-',linewidth=2, markeredgecolor='r',markeredgewidth=1, markersize=6,label='Free Lys (simulated)')
plt.plot(tspan_medium,A1_P_ratio_simu[:,1], 'b-',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=6,label='Protein (simulated)')
plt.plot(tspan_medium,A1_P_ratio_simu[:,2:-1], 'b-',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=6)
plt.plot(tspan_medium,A1_P_ratio_simu[:,len_P_data,], 'g-',linewidth=2, markeredgecolor='g',markeredgewidth=1, markersize=6,label='UIP')
plt.plot(tspan, A1_ratio,     'rs',linewidth=2, markeredgecolor='r',markeredgewidth=2, markersize=8,label='Free Lys(data)')
plt.plot(justify(tspan_data[:,1], invalid_val=-1, axis=0, side='left'),justify(P_ratio[:,1],invalid_val=np.nan, axis=0, side='left'), 'bo',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=8,label='Protein(data)')
plt.plot(justify(tspan_data, invalid_val=-1, axis=0, side='left'),justify(P_ratio,invalid_val=np.nan, axis=0, side='left'), 'bo',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=5)
plt.plot(tspan_medium,A1_P_ratio_simu[:,0], 'r-',linewidth=2, markeredgecolor='r',markeredgewidth=1, markersize=6)
plt.title('Original data and simulations',FontSize = 16)
plt.xlabel('Time (days)', FontSize = 14); plt.ylabel('Fraction of  \'light\'  Lys ', FontSize=14)
plt.legend(loc='best'); plt.ylim(0, 1.03);# plt.show()
plt.savefig('Results/simulated_data.jpeg', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', papertype='None',transparent='False', pad_inches=0.1, frameon='None', metadata='None')

