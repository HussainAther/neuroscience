import numpy as np
import matplotlib.pyplot as plt
from PIL import Image,ImageEnhance
import os
from __future__ import print_function
import argparse
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torchvision import datasets, transforms
import time


"""
Convolutional neural netowkr using torch.
"""

imageDim = 28 # image dimension
numClasses = 10 # number of classes
filterDim = 9 # number of dimensions of filter
numFilters = 20
poolDim = 2

