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

import cal_statics

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

class lorenz96_score:
  def __init__(self):
    print('Call ..... Constract of lorenz96_score')

  def csv2list(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")
  
  def accuracy_rmse_func(self, true_list:list, asses_list:list) -> list:
    rmse = cal_statics.rmse(true_list, asses_list)
    return rmse

  def rmse_draw(self, timestep:int, rmse_da:list) -> None:
    fig = plt.figure()
    ax1 = fig.subplots()
    
    time_list = np.arange(timestep)
    sns.set_style('whitegrid')
    ax1.plot(time_list, rmse_da, ls="--", color='r', label='Data assim.')
    ax1.scatter(time_list, rmse_da, marker='o', color='g', alpha=0.5)
    #ax1.scatter(self.obs_time_list, rmse_obs, marker='o', color='g', s=20, alpha=0.5, edgecolor='k', label='Obs.')

    ax1.set_xlabel('day')
    ax1.set_ylabel('RMSE')
    ax1.set_ylim(0, 10)
    ax1.set_xlim(0, 40)
    ax1.set_title('Lorenz(1996) RMSE, Data assim. method: ' + da_method, loc='left')

    plt.grid()
    plt.legend()
    plt.show()

  def draw_hovmoller(self, nx:int, timestep:int, lorenz_data:np.ndarray) -> None:
    fig, ax = plt.subplots()
    xs = np.arange(1,nx+1)
    ts = np.arange(timestep)
    levels  = np.arange(-6.0, 6.0, 0.5) 
    cmap = plt.get_cmap('bwr')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cmap = plt.contourf(xs, ts, lorenz_data, levels, cmap=cmap, norm=norm, extend='both')
    ax.set_ylim(ax.get_ylim()[::-1])
    fig.colorbar(cmap, ax=ax)
    ax.set_title('contourf with levels', loc='left')

    plt.show()


if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  #---------------------------------------------------------- 
  nx = 24
  stepday = 40
  timestep = stepday*4+1
  
  # OBS
  obs_xintv, obs_tintv = 2, 2
  ny = int(nx/obs_xintv)
  obs_timestep = int(timestep/obs_tintv) + 1

  outdir = './output/lorenz96'
  da_method = 'KF'

  data_true_path = outdir + '/normal_true_score.csv'
  data_NoDA_path = outdir + '/normal_NoDA_score.csv'
  data_DA_path   = outdir + '/normal_'+ da_method + '_DA_score.csv'
  data_obs_path  = outdir + '/normal_obs_score.csv'

  #---------------------------------------------------------- 
  # +++ class set
  # > lorenz96 cal. score
  #---------------------------------------------------------- 
  lz = lorenz96_score()

  #---------------------------------------------------------- 
  # +++ reading func.
  #---------------------------------------------------------- 
  true_score = lz.csv2list(data_true_path).reshape(timestep, nx)
  noda_score = lz.csv2list(data_NoDA_path).reshape(timestep, nx)
  da_score   = lz.csv2list(data_DA_path).reshape(timestep, nx)
  obs_score  = lz.csv2list(data_obs_path).reshape(obs_timestep, ny)
  
  #lz.draw_hovmoller(nx, timestep, true_score)
  lz.draw_hovmoller(nx, timestep, noda_score)
  lz.draw_hovmoller(nx, timestep, da_score)

  #---------------------------------------------------------- 
  # +++ RMSE func.
  #---------------------------------------------------------- 
  rmse_da = []
  for it in range(timestep):
    rmse = lz.accuracy_rmse_func(true_score[it], da_score[it])
    rmse_da.append(rmse)

  lz.rmse_draw(timestep, rmse_da)


  

  