import numpy as np
import matplotlib.pyplot as plt

"""
The Hodgkin-Huxley (HH Hodking Huxley hodgking huxley) model was developed on the basis of 
data from the giant squid axon. The Connor-Stevens model provides an alternate description 
of action-potential general using fast Na+, delayed rectifier K+, and leakage 
conductances (as the HH model does), but the Na+ and K+ conductances have properties 
somewhat different. The C-S model also has an extra K+ conductance (known as the A-current) that is transient.
"""

dt = 0.000005 # step size for time
tmax = 2 # max time

iclamp_flag = 1 # if this is 1, run under current clamp conditions
vclamp_flag = 0 # otherwise this should be 1, for voltage clamp conditions

istart = 0.5 # time applied current starts
ilength = 1 # length of applied current pulse
Ie = 2.405e-7 # magnitude of applied current pulse
Ie = 0.78e-7 # threshold for constant spiking with no A-current

vstart = 0.25 # time to step voltage
vlength = 0.5 # length of voltage step
V0 = -0.080 # initial voltage before step
Ve = -0.040 # value of votage stepped to

V_L = -0.070 # leak reversal potential
E_Na = 0.055 # reversal for sodium channels
E_K = -0.072 # reversal for potassium channels
E_A = -0.075 # reversal for A-type potassium channels

g_L = 3e-6 # specific leak conductance
g_Na = 1.2e-3 # specific sodium conductance
g_K = 2e-4 # specific potassium conductance
g_A = 4.77e-4 # specific A-tpe potassium conductance
g_A = 0.0 # if g_A is zero it switches off the A-current
cm = 10e-9 # specific membrane capacitance

t = range(0, tmax, dt) # time vector
V = np.zeros(len(t)) # voltage vector

if iclamp_flag: # i.e. if in current-clamp mode
    V[0] = V_L # set the inititial value of voltage

n = np.zeros(len(t)) # n: potassium activation gating variable
m = np.zeros(len(t)) #  sodium activation gating variable
h = np.zeros(len(t)) #  h: sodim inactivation gating variable

a = np.zeros(len(t)) # A-current activation gating variable
b = np.zeros(len(t)) # A-current inactivation gating variable

Iapp = np.zeros(len(t)) # Applied current, relevant in current-clamp mode
if iclamp_flag: # i.e. if in current-clamp mode
    for i in range(round(istart/dt), round((istart+ilength)/dt)+1) # make non-zero for duration of current pulse
        Iapp[i-1] = Ie

Vapp = np.zeros(len(t)) # Applied voltage, relevant in voltage-clamp mode
if vclamp_flag:
    for i in range(round(vstart/dt)) # make V0 before pulse
        Vapp[i-1] = V0
    for i in range(round(vstart/dt), round((vstart+vlength)/dt)+1) # make Ve for duration of voltage pulse
        Vapp[i-1] = Ve
    for i in range(round((vstart+vlength)/dt), length(Vapp)+1) # make V0 following pulse
        Vapp[i-1] = V0

# total current
Itot = np.zeros(len(t)) # in case we want to plot and look at the total current

# now see how things change over time
for i in range(2, len(t)+1):
    I_L = g_L * (V_L-V[i-2])
    
    Vm = V[i-2] * 1000 # converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
    
    """ 
    Sodium and potassium gating variables are defined by the
    voltage-dependent transition rates between states, labeled alpha and
    beta. Written out from Dayan/Abbott, units are 1/ms.
    """
    
    if Vm == -29.7:
        alpha_m = 0.38/0.1
    else
        alpha_m = 0.38 * (Vm + 29.7)/(1 - exp(-0.1 * (Vm + 29.7)))
    beta_m = 15.2 * exp(-0.0556 * (Vm + 54.7))

    alpha_h = 0.266 * exp(-0.05 * (Vm + 48))
    beta_h = 3.8/(1 + exp(-0.1 * (Vm + 18)))

    if Vm == -45.7:
       alpha_n = 0.02/0.1
    else
        alpha_n = 0.02 * (Vm + 45.7)/(1 - exp(-0.1 * (Vm + 45.7)))
    beta_n = 0.25 * exp(-0.0125 * (Vm + 55.7))
 
    """
    From the alpha and beta for each gating variable we find the steady
    state values (_inf) and the time constants (tau_) for each m,h and n.
    """

    tau_m = 1e-3 / (alpha_m + beta_m) # time constant converted from ms to sec
    m_inf = alpha_m / (alpha_m + beta_m)

    tau_h = 1e-3/(alpha_h + beta_h) # time constant converted from ms to sec
    h_inf = alpha_h/(alpha_h + beta_h)

    tau_n = 1e-3/(alpha_n + beta_n) # time constant converted from ms to sec
    n_inf = alpha_n/(alpha_n + beta_n)

    m[i-1] = m[i-2] + (m_inf-m[i-2]) * dt/tau_m # Update m

    h[i-1] = h[i-2] + (h_inf-h[i-2]) * dt/tau_h # Update h

    n[i-1] = n[i-2] + (n_inf-n[i-2]) * dt/tau_n #Update n
    
    """
    For the A-type current gating variables, instead of using alpha and
    beta, we just use the steady-state values a_inf and b_inf along with
    the time constants tau_a and tau_b that are found empirically
    (Dayan-Abbott, p. 224)
    """

    a_inf = (0.0761 * exp(0.0314 * (Vm + 94.22))/(1 + exp(0.0346 * (Vm + 1.17))))^(1/3.0)
    tau_a = 0.3632*1e-3 + 1.158e-3/(1 + exp(0.0497 * (Vm + 55.96)))

    b_inf = (1/(1 + exp(0.0688 * (Vm + 53.3))))^4
    tau_b = 1.24e-3 + 2.678e-3/(1 + exp(0.0624 * (Vm + 50)))

    a[i-1] = a[i-2] + (a_inf-a[i-2])*dt/tau_a # Update a
    b[i-1] = b[i-2] + (b_inf-b[i-2])*dt/tau_b # Update b

    I_Na = g_Na*m[i-1]*m[i-1]*m[i-1]*h[i-1]*(E_Na-V[i-2]) # total sodium current

    I_K = g_K*n[i-1]*n[i-1]*n[i-1]*n[i-1]*(E_K-V[i-2]) # total potassium current

    I_A = g_A*a[i-1]*a[i-1]*a[i-1]*b[i-1]*(E_A-V[i-2]) # total A-type current

    Itot[i-2] = I_L+I_Na+I_K+I_A+Iapp[i-2] # total current is sum of leak + active channels + applied current

    V[i-1] = V[i-2] + Itot[i-2]*dt/cm # Update the membrane potential, V.

    if vclamp_flag: # if we are using voltage clamp
        V[i-1] = Vapp[i-1] # ignore the voltage integration and set V to be the applied voltage

plt(t,V)
plt.show()
