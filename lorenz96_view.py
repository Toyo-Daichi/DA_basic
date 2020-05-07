#-*- coding: utf-8 -*-
"""
Created on 2020.5.7
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './lib'))

import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class lorenz96_score:
  def __init__(self, path:str, nx:int, nt:int):
    self.lorenz_data = self.csv2list(path).reshape(nt, nx)
    self.xs = np.arange(1, nx+1, step=1)
    self.time_list = np.arange(nt, step=1)
    print('Call ..... Constract of lorenz96_score')

  def csv2list(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")
  
  def lorenz96_hovmoller(self) -> None:
    ax = plt.gca()
    cmap = plt.contourf(self.xs, self.time_list, self.lorenz_data, cmap=cm.coolwarm, extend='both')
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.show()

if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  nx = 24
  timestep = 40*4+1
  
  outdir    = './output/'
  da_method = ''
  data_path = outdir + 'lorenz96' + da_method + '.csv'

  #---------------------------------------------------------- 
  # +++ class set
  # > lorenz96 cal. score
  score = lorenz96_score(data_path, nx, timestep)
  print(score.lorenz_data.shape)
  print(score.xs.shape)
  print(score.time_list.shape)
  #---------------------------------------------------------- 
  # +++ draw func.
  # > hovmoller diagram
  score.lorenz96_hovmoller()