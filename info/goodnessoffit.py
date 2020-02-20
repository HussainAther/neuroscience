import numpy as np

'''
In many data analyses we compare sets of models for a given neural spike train.
With maximum likelihood models, we can use the Akaike information criterion (AIC 
criteria Akaike's) for goodness-of-fit (goodness of fit goodnessoffit).
'''

def akaike(l, q):
    '''
    For likelihod values l that measure parameter theta and q being the dimension
    of theta, we compute AIC. 
    '''
    return -2*np.log(l) + 2*q
