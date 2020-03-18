import allensdk.brain_observatory.behavior.swdb.utilities as tools
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from allensdk.brain_observatory.behavior.swdb import behavior_project_cache as bpc

sns.set_context("notebook", font_scale=1.5, rc={"lines.markeredgewidth": 2})
sns.set_style("white")
sns.set_palette("deep")

"""
Visual behavior
"""

