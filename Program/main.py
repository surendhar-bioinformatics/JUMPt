import numpy as np, scipy.optimize as so, scipy.integrate as si, pandas as pd, de, numpy.linalg as nl, coupledOde as co, scipy.sparse as ss, generateUnobs as gu, matplotlib.pyplot as plt, time, subprocess, os
from datetime import datetime 
from pathlib import Path
import tkinter


def runMatlab( scriptFile ):
    """
    Run matlab on scriptFile and return when finished.
    """
    retval = subprocess.run( ['matlab', '-wait', '-nojvm','-nodisplay', '-nosplash',  '-r', scriptFile ] )
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


################################################################### 
#================== 1.Load the original data ======================
################################################################### 
# Ask the user to select a single file name.
#Inputfile   = "test_data.xlsx" # Assign spreadsheet filename to `file`
Inputfile = input('Enter your input data file  (if it doesnt exist in current folder, enter with complete path;  test data file name is "test_data.xlsx"):') 
#root = tkinter.Tk()
#my_filetypes = [('all files', '.xlsx'), ('Excel file', '.txt')] # Build a list of tuples for each file type the file dialog should display
#Inputfile = tkinter.filedialog.askopenfilename(initialdir=os.getcwd(), title="Please select a file:", filetypes=my_filetypes)   
#root.destroy()
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
plt.plot(tspan, A1_ratio,     'r-s',linewidth=2, markeredgecolor='r',markeredgewidth=2, markersize=10,label='Free Lys')
plt.plot(justify(tspan_data[:,1], invalid_val=-1, axis=0, side='left'),justify(P_ratio[:,1],invalid_val=np.nan, axis=0, side='left'), 'b-o',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=10,label='Protein')
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
print('\n\n Completed indiviual protein fitting and results were saved')

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

res_UIP = so.minimize_scalar( noFun, bounds=[0,20], 
                          method='bounded', options={'disp':3,'maxiter':2000} )
gamma_UIP = res_UIP.x
Theta_UIP = gu.generateSys( t, df[:,0], gamma_UIP )
dTheta_UIP= gamma_UIP*(f[:,0] - Theta_UIP)
c        = np.dot(ut.T,-dTheta_UIP)/(nl.norm(dTheta_UIP)**2)*A1
EtaP_forDeriFit = np.hstack((EtaP, c))

gamma     = np.hstack((gammaUt[list(range(0,len_P_data))].reshape((-1,)), gamma_UIP, gammaUt[len_P_data])); #gamma =abs(gamma)
integrator = co.coupledOde( len(EtaP_forDeriFit), EtaP_forDeriFit, A1,
                            np.ones(len(EtaP_forDeriFit)), np.ones(1), np.ones(1)*.05, gamma )


refft   = np.hstack((A1_ratio, P_ratio[:,:])); 
refft   = np.vstack([x.reshape((1,-1)) for x in refft])
rv      = integrator.integrate(tspan)
idx     = list(range(len(EtaP_forDeriFit)-1)) + [-1]                 
residual= rv[idx,:] - refft[:,list(range(1,len(EtaP_forDeriFit)))+[0]].T
mask    = (1 - np.isnan(residual)).nonzero()
flatResidual = residual[mask]; ResError = nl.norm(flatResidual)**2;
ResError_temp = .001

while ResError_temp <= ResError:
    print('\n Working to minize the residual error ; Thanks for your patience ...... \n')
    res     = so.minimize( lambda x: nl.norm(np.dot(G,x) - MTb.T,2)**2 + nl.norm(x*(x < 0),2)**2, np.array(np.abs(gammaUt_temp)).reshape((-1,)), jac=lambda x: (np.dot(G,x) - MTb.T) + (x*(x < 0)), method='tnc', options={'maxiter':126} )
    res     = so.minimize( lambda x: nl.norm(np.dot(G,x) - MTb.T,2)**2, res.x, jac=lambda x: (np.dot(G,x) - MTb.T), method='trust-constr', options={'maxiter':126}, bounds=[(0,None)]*G.shape[0] )
    gammaUt_temp = res.x.reshape((-1,1))
    ut      = np.array(gammaUt_temp[A.shape[1]:]).reshape((-1,1))
    ut[0] = 0
    
    def noFun( gamma ):
        fP = gu.generateSys( t, df[:,0], gamma )
        return np.sin(np.arccos((np.dot(gamma*(f[:,0] - fP),ut)/(nl.norm(gamma*(f[:,0] - fP))*nl.norm(ut)))))
    
    res_UIP = so.minimize_scalar( noFun, bounds=[0,20], 
                              method='bounded', options={'maxiter':2000} )
    gamma_UIP_temp  = res_UIP.x
    Theta_UIP       = gu.generateSys( t, df[:,0], gamma_UIP_temp )
    dTheta_UIP      = gamma_UIP_temp*(f[:,0] - Theta_UIP)
    c               = np.dot(ut.T,-dTheta_UIP)/(nl.norm(dTheta_UIP)**2)*A1
    EtaP_forDeriFit_temp = np.hstack((EtaP, c))
    gamma_temp     = np.hstack((gammaUt[list(range(0,len_P_data))].reshape((-1,)), gamma_UIP_temp, gammaUt[len_P_data])); #gamma =abs(gamma)

    
    integrator = co.coupledOde( len(EtaP_forDeriFit_temp), EtaP_forDeriFit_temp, A1,
                                np.ones(len(EtaP_forDeriFit_temp)), np.ones(1), np.ones(1)*.05, gamma_temp )
    
    rv      = integrator.integrate(tspan)
    idx     = list(range(len(EtaP_forDeriFit_temp)-1)) + [-1]                 
    residual= rv[idx,:] - refft[:,list(range(1,len(EtaP_forDeriFit_temp)))+[0]].T
    mask    = (1 - np.isnan(residual)).nonzero()
    flatResidual = residual[mask]; ResError_temp = nl.norm(flatResidual)**2; 
    print('\n Cuurent Residual Error = ', ResError_temp,' Previous Residual Error = ', ResError)
    if ResError_temp < ResError:
        ResError    = ResError_temp;
        gammaUt     = gammaUt_temp
        gamma_UIP   = gamma_UIP_temp
        gamma       = gamma_temp
        EtaP_forDeriFit = EtaP_forDeriFit_temp
        


integrator = co.coupledOde( len(EtaP_forDeriFit), EtaP_forDeriFit, A1,
                            np.ones(len(EtaP_forDeriFit)), np.ones(1), np.ones(1)*.05, gamma )

rv      = integrator.integrate(tspan)
idx     = list(range(len(EtaP_forDeriFit)-1))+[-1]                  
residual= rv[idx,:] - refft[:,list(range(1,len(EtaP_forDeriFit)))+[0]].T
mask    = (1 - np.isnan(residual)).nonzero()
flatResidual = residual[mask]; ResError = nl.norm(flatResidual)**2;
print('\n \n Residual Error (LSE) after Global fitting of Proteins and Lys = ', ResError)

rv_medium      = integrator.integrate(tspan_medium)
idx3 = [len_P_data+1] + list(range(0,len_P_data+1)); 
rv_medium  = rv_medium[idx3, :]
half_life = np.log(2)/gamma; gamma3 = (np.vstack((gamma,half_life))).T
row     = [P_data.iloc[i,0] for i in range(1, len_P_data+2)]; row = row +['UIP'];
Gamma_DeriFit = pd.DataFrame(gamma3,index = row, columns = ['degredation_rate', 'HalfLife']); 
writer  = pd.ExcelWriter('Results/Glob_fit_res.xlsx', engine='xlsxwriter'); # Create a Pandas Excel writer using XlsxWriter as the engine.
Gamma_DeriFit.to_excel(writer, sheet_name='GlobParams') # Convert the dataframe to an XlsxWriter Excel object.

traj_Lys_P_DeriFit = pd.DataFrame(rv_medium.T,columns = row, ); 
traj_Lys_P_DeriFit.to_excel(writer, sheet_name='traj_Lys_P',  index=False) # Convert the dataframe to an XlsxWriter Excel object.
writer.save()# Close the Pandas Excel writer and output the Excel file.
writer.close()

Matlab_GlobOpti = input('Do you wish to run Global optimization in MATLAB ? (Yes OR No):') 
################################################################### 
###============== Global optimization for degrdation rates======##
################################################################### 
if Matlab_GlobOpti == "Yes" or  Matlab_GlobOpti == "yes" or  Matlab_GlobOpti == "y" or Matlab_GlobOpti == "Y" :
    print('\n Calling MATLAB program to fit the Lys and all proteins Globally ...')
    runMatlab( 'PT_Glob_Fitting' )
    print('\n Completed Global fitting and results were saved')

#########################################################################################################################  
###============ Generating figure with Global optimized parameters  and saving the tragectries of all proteins  ======##
#########################################################################################################################  

Inputfile   = Path("Results\Glob_fit_res.xlsx") # Assign spreadsheet filename to `file`
A1_P_ratio_simu  = pd.read_excel(Inputfile,sheet_name= 'traj_Lys_P') # Load spreadsheet
A1_P_ratio_simu  = A1_P_ratio_simu.values
plt.figure(2)
plt.plot(tspan_medium,A1_P_ratio_simu[:,0], 'r-',linewidth=2, markeredgecolor='r',markeredgewidth=1, markersize=6,label='Free Lys (simulated)')
plt.plot(tspan_medium,A1_P_ratio_simu[:,1], 'b-',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=6,label='Protein (simulated)')
plt.plot(tspan_medium,A1_P_ratio_simu[:,2:-1], 'b-',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=6)
plt.plot(tspan_medium,A1_P_ratio_simu[:,len_P_data+1], 'c-',linewidth=2, markeredgecolor='c',markeredgewidth=1, markersize=6,label='UIP')
plt.plot(tspan, A1_ratio,     'rs',linewidth=2, markeredgecolor='r',markeredgewidth=2, markersize=8,label='Free Lys(data)')
plt.plot(justify(tspan_data[:,1], invalid_val=-1, axis=0, side='left'),justify(P_ratio[:,1],invalid_val=np.nan, axis=0, side='left'), 'bo',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=8,label='Protein(data)')
plt.plot(justify(tspan_data, invalid_val=-1, axis=0, side='left'),justify(P_ratio,invalid_val=np.nan, axis=0, side='left'), 'bo',linewidth=2, markeredgecolor='b',markeredgewidth=1, markersize=5)
plt.plot(tspan_medium,A1_P_ratio_simu[:,0], 'r-',linewidth=2, markeredgecolor='r',markeredgewidth=1, markersize=6)
plt.title('Original data and simulations',FontSize = 16)
plt.xlabel('Time (days)', FontSize = 14); plt.ylabel('Fraction of  \'light\'  Lys ', FontSize=14)
plt.legend(loc='best'); plt.ylim(0, 1.03);# plt.show()
plt.savefig('Results/simulated_data.jpeg', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', papertype='None',transparent='False', pad_inches=0.1, frameon='None', metadata='None')

