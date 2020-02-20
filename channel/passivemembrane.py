import matplotlib.pyplot as plt
import numpy as np

'''
Passive membrane subject to current pulse solved with
Backward Euler method.
'''

def bepswI(dt, Tfin):
    '''
    Bakcwards Euler passive membrane 
    '''
    VCl = -68 # mV source chloride voltage at beginning
    A = 4*np.pi*1e-6 # cm^2 patch area
    Cm = 1 # micro F/cm^2 capacitance per area
    gCl = .3 # mS/cm^2 conductance
    tau = Cm/gCl # ms time constant
    Nt = round(1+Tfin/dt) # time steps
    v = np.zeros(Nt) # voltage
    t = np.zeros(Nt) # time
    ICl = np.zeros(Nt) # chloride current
    IC = np.zeros(Nt) # capactive current
    v[0] = VCl
    for j in range(1, Nt):
        t[j] = j*dt
        Istim = (t[j]>2)*(t[j]<22)*1e-5 # 10 pA 20 ms pulse
        v[j] = (v[j-1] + dt*(VCl/tau + Istim/A/Cm))/(1+dt/tau) # backward euler
        ICl[j] = gCl*(v[j]-VCl) 
        IC[j] = Istim/A - ICl[j] 
    plt.subplot(1,2,1)
    plt.plot(t, v)
    plt.xlabel('t  (ms)',fontsize=14)
    plt.ylabel('V  (mV)',fontsize=14)
    plt.subplot(1,2,2)
    plot(t, IC, 'r', t, ICl, 'k')
    plt.legend(['$I_C$','$I_{Cl}$'], loc='best')
    plt.xlabel('t  (ms)',fontsize=14)
    plt.ylabel('I  $(\mu A/cm^2)$',fontsize=14)
    plt.show()
