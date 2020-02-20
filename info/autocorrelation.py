import numpy as np

'''
We use the autocorrelation of a time series function to look for repeating patterns
or serial correlation (correlation between signal at a given time and a later time).
'''

def autocorr(x):
    '''
    Numpy's correlate a function with itself.
    '''
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]
