"""
Euclidian spike classifier. To quantify how well a spike train discriminates 
stimulus patterns, one first needs to estimate the optimal stimulus feature 
for eliciting spikes. Here, we describe the application of a Euclidian pattern 
classifier to this problem. First, the spike train x and stimulus s are binned to
give us a max of one spike per bin. r_t takes the value of 1 if the time bin ending
at t has a spike and the value 0 if it doesn't. Segments s_t ending at time t and
comprising of ~100 bins prior to time t are assigned to one of two ensembles 
P(s|r=1) and P(s|r=0) depending on whether there was a spike or not. Compute feature f
from the means m1 and m= of the conditional distributions such taht f = m1 - m0.
"""

def epc(bins):
    """
    For array bins of 0s and 1s corresponding to spike trains, compute the Euclidian
    pattern classifier.
    """
    r_t = [] # rate array
    one = False # whether we are viewing spikes in the iteration
    currscore = 0 # current score
    for i in bins:
        if i == 1:  
            currscore += 1
            if one == False:
                one = True
                r_t.append(currscore)
                currscore = 0 
            else:
                currscore += 1  
        else:
            r_t.append(0)
    zeroindices = [i for i, x in enumerate(bins) if x == 0]
    oneindices = [i for i, x in enumerate(bins) if x == 1]
    m0 = len(zeroindices)/len(bins)
    m1 = len(oneindices)/len(bins)
    return r_t, m0, m1
