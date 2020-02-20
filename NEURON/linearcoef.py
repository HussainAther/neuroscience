import matplotlib.pyplot as plt
import numpy as np

'''
Find linear coefficients of conductance for a 
potassium channel.
'''

gkbar_range = numpy.arange(.001, 0.012, 0.001)
max_voltages = []
for gkbar in gkbar_range:
    soma.gkbar_hh = gkbar
    neuron.h.run()
    max_voltages.append(max(voltage))

linear_coef = np.polyfit(gkbar_range, max_voltages, 1)
print('Linear equation max_voltage = %f*gkbar + %f' % tuple([x for x in linear_coef]))
plt.plot(gkbar_range, max_voltages)
plt.xlabel('gkbar (S/cm2)')
plt.ylabel('Maximum AP voltage')
for xs in [0.1, 0.15]:
    plt.axvline(x=xs, color='r')
plt.show()

