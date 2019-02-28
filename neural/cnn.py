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
numFilters = 20 # number of filters
poolDim = 2  # number of pools uesd in the final step

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.layer1 = nn.Sequential(
            nn.Conv2d(1,numFilters, kernel_size=filterDim, stride=1),
            nn.Sigmoid(),
            nn.MaxPool2d(kernel_size=poolDim, stride=2))
        #28x28x1 => (28-filterDim+1)/2
        dim=(imageDim-filterDim+1)/poolDim
        
        self.fc1 = nn.Linear(int(dim*dim*numFilters), 10)
    
    
    
    def forward(self, x):
        #x = F.relu(self.affine1(x))
        #action_scores = self.affine2(x)
        out = self.layer1(x)
        out = out.reshape(out.size(0), -1)
        out = self.fc1(out)
        #out = self.BatchNorm2(out)
        #action_scores = self.fc2(out)
        return F.softmax(out, dim=1)


def train(args, model, device, train_loader, optimizer, epoch):
    model.train()
    for batch_idx, (data, target) in enumerate(train_loader):
        data, target = data.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = F.nll_loss(output, target)
        loss.backward()
        optimizer.step()
        if batch_idx % args.log_interval == 0:
            print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                epoch, batch_idx * len(data), len(train_loader.dataset),
                100. * batch_idx / len(train_loader), loss.item()))


