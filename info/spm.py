import numpy as np
import scipy.sparse as sps

"""
Statistical Parametric Mapping refers to the construction and assessment of 
spatially extended statistical processes used to test hypotheses about functional 
imaging data. These ideas have been instantiated in software that is called SPM.
"""

def spm_vec(X):
    """
    SPM Statistical parametric mapping. Vectorize a numeric, cell or structure array X.
    """
    for v_it, v in enumerate(X):
        if v_it == 0:
            try:
                nr,nc = v.shape
                vX = np.reshape(v,[nr*nc,1])
            except:
                vX = np.array([v][:,np.newaxis])
        else:
            try:
                nr,nc = v.shape
                vX = np.concatenate([vX,np.reshape(v,[nr*nc,1])])
            except:
                vX = np.concatenate([vX,np.array([v])[:,np.newaxis]])            
    return np.squeeze(vX)[:,np.newaxis]


def spm_unvec(Xflat, X):
    """
    SPM Statistical parametric mapping from flat map Xflat and given X values X onto a list.
    """
    count = 0
    vXlist = []
    for v_it,v in enumerate(X):
        try:
            nr,nc = v.shape
            rng = np.arange(cnt,nr*nc)
            vX = np.reshape(Xflat[rng],[nr,nc])
            vXlist.append(vX)
            count+=1
        except: 
            vXlist.append(v)
    return vXlist

def erp(x, u, P, M, returnJ=True, returnD=True):
    """
    Event-related potential for state vector x with the following:
    x[1] - voltage (spiny stellate cells)
    x[2] - voltage (pyramidal cells) +ve
    x[3] - voltage (pyramidal cells) -ve
    x[4] - current (spiny stellate cells)    depolarizing
    x[5] - current (pyramidal cells)         depolarizing
    x[6] - current (pyramidal cells)         hyperpolarizing
    x[7] - voltage (inhibitory interneurons)
    x[8] - current (inhibitory interneurons) depolarizing
    x[9] - voltage (pyramidal cells)
    with P as true connectivity parameters, u as potential energy, and M as the dictionary of
    types of models. returnJ and returnD are for returning Jacobian and delays, respectively.
    Return f (dx(t)/dt = f(x(t))), J (df(t)/dx(t)), and D (delay operator dx(t)/dt)
    """
    n = len(P["A"][0]) # number of sources
    x = spm_unvec(x, M["x"]) # extract neuronal states

    E = np.array([1., 1/2., 1/8.,])*32 # extrinsic rates (forward, backward, lateral)
    G = np.array([1, 4/5., 1/4., 1/4*128.]) # intrinsic rates (g1 g2 g3 g4)    
    D = np.array([2, 16])  #  propogation delays (intrinsic, extrinsic)
    H = np.array([4, 32])  #  receptor densities (excitatory, inhibitory)
    T = np.array([8, 16])  #  synaptic constants (excitatory, inhibitory)
    R = np.array([2, 1])/3  # parameters of static nonlinearity

    try:  # test for free parameters on intrinsic connections
        G = G*np.exp(P["H"])
    except: _
        G = np.ones([n,1])*G
    
    if "pF" in M: # which parameters are fixed
        try: E = M["pF"]["E"] 
        except: _
        try: G = M["pF"]["H"] 
        except: _
        try: D = M["pF"]["D"] 
        except: _
        try: H = M["pF"]["G"] 
        except: _
        try: T = M["pF"]["T"] 
        except: _  
        try: R = M["pF"]["R"] 
        except: _  

    A = [] # exponential transofrm to ensure positivity constraints
    A[0] = np.exp(P["A"][0])*E[0]
    A[1] = np.exp(P["A"][1])*E[1]    
    A[2] = np.exp(P["A"][2])*E[2]
    C = np.exp(P["C"])

    # Intrinsic connectivity and parameters
    Te = T[0]/1000*np.exp(P["T"][0]) # excitatory time constants
    Ti = T[1]/1000*np.exp(P["T"][1]) # inhibitory time constants
    He = H[0]*np.exp(P["G"][0]) # excitatory receptor density
    Hi = H[1]*np.exp(P["G"][1]) # inhibitory receptor density

    R = R*np.exp(P["S"]) # pre-synaptic inputs
    S = 1./(1 + np.exp(-R[0]*(x[0] - R[1]))) - 1./(1 + np.exp(R[0]*R[1])) 
  
    if "u" in M: # endogenous input
        U = u[:]*64
    else: # exogenous input
        U = C*u[:]*2
   
    f = [] # output velocity
    # Supragranular layer (inhibitory interneurons) with voltage and depolarizing current
    f[6] = x[7]
    f[7] = (He*( (A[1] + A[2]) * S[8] + G[2]*S[8]) - 2*x[7] - x[6]/Te)/Te

    # Granular layer (spiny stellate cells) with voltage and depolarizing current
    f[0] = x[3]
    f[3] = (He*((A[0] + A[2])*S[8] + G[0]*S[8] + U) - 2*x[3] - x[0]/Te)/Te

    # Infra-granular layer (pyramidal cells) with depolarizing current
    f[1] = x[4]
    f[4] = (He*((A[1] + A[2])*S[8] + G[1]*S[0]) - 2*x[4] - x[1]/Te)/Te

    # Infra-granular layer (pyramidalcells) with hyperpolarizing current
    f[2] = x[5]
    f[5] = (Hi*G[3]*S[6] - 2*x[5] - x[2]/Ti)/Ti

    # Infra-granular layer (pyramidal cells) with voltage
    f[8] = x[4] - x[5]

    f = spm_vec(f) # vectorize 

    if returnJ == False and returnD == False:
        return f
    J = np.gradient([M["f"], x, u, P, M, 1]) # Jacobian gradient
    
    # Delays
    De = D[1] * np.exp(P["D"])/1000
    Di = D[0] / 1000
    De = np.array*([1-sps.eye(n,n)]) *De # sps.eye returns sparse matrix with ones in diagnol
    Di = np.array*([1-sps.eye(9,9)]) *Di
    De = np.kron(np.ones([9,9]),De) # Kronecker product
    Di = np.kron(Di,sps.eye([n,n]))
    D = Di + De
    
    # Q matrix from Jacobian
    Q = np.linalg.inv(sps.eye(len(J)) + D*J)
    return f, j, D
