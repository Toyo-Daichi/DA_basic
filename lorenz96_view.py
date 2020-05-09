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
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

class lorenz96_score:
  def __init__(self, path:str, nx:int, nt:int):
    self.lorenz_data = self.csv2list(path).reshape(nt, nx)
    self.xs = np.arange(1, nx+1, step=1)
    self.ts = np.arange(nt, step=1)
    print('Call ..... Constract of lorenz96_score')

  def csv2list(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")
  
  def lorenz96_hovmoller(self) -> None:
    fig, ax = plt.subplots()
    levels = [ 
                     -1.8, -1.6, -1.4, -1.2, -1.0,
                      1.2,  1.4,  1.6,  1.8,  2.0,
                      2.2,  2.4,  2.6,  2.8,  3.0 ]
    cmap = plt.get_cmap('bwr')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cmap = plt.contourf(self.xs, self.ts, self.lorenz_data, levels, cmap=cmap, norm=norm, extend='both')
    ax.set_ylim(ax.get_ylim()[::-1])
    fig.colorbar(cmap, ax=ax)
    ax.set_title('contourf with levels', loc='left')

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

  #---------------------------------------------------------- 
  # +++ draw func.
  # > hovmoller diagram
  score.lorenz96_hovmoller()