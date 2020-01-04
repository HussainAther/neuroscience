import neurodynex.ojas_rule.oja as oja

"""
Oja's learning rule
"""

# Generate data points.
cloud = oja.make_cloud()

# Learn weights and return timecourse. 
wcourse = oja.learn(cloud)  
