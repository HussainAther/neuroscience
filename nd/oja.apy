import neurodynex.ojas_rule.oja as oja
import matplotlib.pyplot as plt

"""
Oja's learning rule
"""

# Generate data points.
cloud = oja.make_cloud()

# Learn weights and return timecourse. 
wcourse = oja.learn(cloud)  

# Plot.
plt.scatter(cloud[:, 0], cloud[:, 1], marker=".", alpha=.2)
plt.plot(wcourse[-1, 0], wcourse[-1, 1], "or", markersize=10)
plt.axis("equal")
plt.figure()
plt.plot(wcourse[:, 0], "g")
plt.plot(wcourse[:, 1], "b")
print("The final weight vector w is: ({},{})".format(wcourse[-1,0],wcourse[-1,1]))
