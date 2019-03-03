import numpy as np

"""
The Hodgkin-Huxley model was developed on the basis of data from the giant squid axon.
The Connor-Stevens model provides an alternate descirption of action-potential general using
fast Na+, delayed rectifier K+, and leakaage conductsances (as the HH model does), but the Na+ and K+
conductances have properties somewhat different. The C-S model als ohas an extra K+ conductance
(know nas the A-current) that is transient.
"""

dt = 0.000005
tmax = 2

iclamp_flag = 1 # if this is 1, run under current clamp conditions
vclamp_flag = 0 # otherwise this should be 1, for voltage clamp conditions

istart = 0.5 # time applied current starts
ilength = 1 # length of applied current pulse
Ie=2.405e-7 # magnitude of applied current pulse
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
    V(1) = V_L # set the inititial value of voltage

n = np.zeros(len(t)) # n: potassium activation gating variable
m = np.zeros(len(t)) #  sodium activation gating variable
h = np.zeros(len(t)) #  h: sodim inactivation gating variable

a = np.zeros(len(t)) # A-current activation gating variable
b = np.zeros(len(t)) # A-current inactivation gating variable

Iapp = np.zeros(len(t)) # Applied current, relevant in current-clamp mode
if iclamp_flag: # i.e. if in current-clamp mode
    for i in range(round(istart/dt), round((istart+ilength)/dt+1)) # make non-zero for duration of current pulse
        Iapp[i] = Ie

Vapp = np.zeros(len(t)) # Applied voltage, relevant in voltage-clamp mode
if ( vclamp_flag )
    for i = 1:round(vstart/dt)          % % make V0 before pulse
        Vapp(i) = V0;
    end
    for i=round(vstart/dt)+1:round((vstart+vlength)/dt) % make Ve for duration of voltage pulse
        Vapp(i) = Ve;
    end
    for i=round((vstart+vlength)/dt):length(Vapp) % make V0 following pulse
        Vapp(i) = V0;
    end
end


def I_m():
    """
    Membrane current for Connor-Stevens Model.
    """

