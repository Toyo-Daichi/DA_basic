#-*- coding: utf-8 -*-
"""
Created on 2020.5.7
@author: Toyo_anlichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './lib'))

import itertools
import numpy as np
import pandas as pd
import seaborn as sns

import cal_statics

from tqdm import tqdm
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

  def rmse_draw(self, timestep:int, rmse_anl:list, rmse_sim:list, rmse_obs:list, obs_tintv:int) -> None:
    fig = plt.figure()
    ax1 = fig.subplots()
    
    time_list = np.arange(timestep)
    obs_time_list = np.arange(obs_tintv, timestep-obs_tintv, obs_tintv)

    sns.set_style('whitegrid')
    ax1.plot(time_list, rmse_anl, ls="--", color='r', label='DA')
    ax1.plot(time_list, rmse_sim, ls="--", label='No DA')
    ax1.scatter(obs_time_list, rmse_obs, marker='*', color='y', s=20, alpha=0.5, edgecolor='k', label='Obs.')

    ax1.set_xlabel('day')
    ax1.set_ylabel('RMSE')
    ax1.set_ylim(0, 10)
    ax1.set_title('Lorenz(1996) RMSE, Data assim. method: ' + da_method, loc='left')

    plt.grid()
    plt.legend()
    plt.show()

  def draw_hovmoller(self, nx:int, timestep:int, lorenz_anlta:np.ndarray) -> None:
    fig, ax = plt.subplots()
    xs = np.arange(1,nx+1)
    ts = np.arange(timestep)
    levels  = np.arange(-8.0, 8.0, 0.5) 
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cmap = plt.contourf(xs, ts, lorenz_anlta, levels, cmap=cmap, norm=norm, extend='both')
    ax.set_ylim(ax.get_ylim()[::-1])
    fig.colorbar(cmap, ax=ax)
    ax.set_title('Lorenz(1996) NX=40, F=8.0', loc='left')

    plt.show()

class lorenz96_errcov:
  """
  一つのデータリストを誤差共分散行列の扱いを行うクラスメソッド
  """
  def __init__(self, path:str):
    self.err_anlta = self.read_error_csv(path)

  def read_error_csv(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")

  def error_heatmap(self, err_anlta:np.ndarray, timestep:int) -> None:
    fig, ax = plt.subplots()
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.set(style='white')
    sns.heatmap(
      data=err_anlta, cmap=cmap, 
      vmax=1.0, vmin=-1.0, center=0,
      square=True, linewidths=.5, 
      cbar_kws={"shrink": .5, "extend": 'both'},
    )

    time = timestep*0.1
    ax.set_title('Back ground Error Cov. - {:.1f} tmp step.'.format(time), loc='center')
    plt.savefig('./figure/error_heatmap_{:.1f}sec.png'.format(time))
    plt.close('all')


if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  #---------------------------------------------------------- 
  nx = 40
  stepday = 40
  timestep = stepday*4+1
  
  # OBS
  obs_xintv, obs_tintv = 1, 1
  ny = int(nx/obs_xintv)
  obs_timestep = int(timestep/obs_tintv)-1
  mems = 50

  outdir = './output/lorenz96'
  da_method = 'EnKF'


  data_true_path = outdir + '/normal_true_score_' + str(nx) + 'n.csv'
  data_sim_path = outdir + '/normal_sim_score_' + str(nx) + 'n.csv'
  data_anl_path   = outdir + '/normal_'+ da_method + str(mems) + 'm_anl_score_' + str(nx) + 'n.csv'
  data_obs_path  = outdir + '/normal_obs_score_' + str(nx) + '.csv'

  data_err_path = outdir + '/Error_matrix_' + da_method + '_' + str(nx) + 'n.csv'

  #---------------------------------------------------------- 
  # +++ class set
  # > lorenz96 cal. score
  #---------------------------------------------------------- 
  lz = lorenz96_score()
  #err = lorenz96_errcov(data_err_path)

  #---------------------------------------------------------- 
  # +++ reading func.
  #---------------------------------------------------------- 
  true_score = lz.csv2list(data_true_path).reshape(timestep, nx)
  noda_score = lz.csv2list(data_sim_path).reshape(timestep, nx)
  anl_score  = lz.csv2list(data_anl_path).reshape(timestep, nx)
  obs_score  = lz.csv2list(data_obs_path).reshape(obs_timestep, ny)

  lz.draw_hovmoller(nx, timestep, true_score)
  lz.draw_hovmoller(nx, timestep, noda_score)
  lz.draw_hovmoller(nx, timestep, anl_score)

  #---------------------------------------------------------- 
  # +++ RMSE func.
  #---------------------------------------------------------- 
  rmse_anl_list = []
  rmse_sim_list = []
  rmse_obs_list = []
  for it in range(timestep):
    rmse_anl = lz.accuracy_rmse_func(true_score[it], anl_score[it])
    rmse_sim = lz.accuracy_rmse_func(true_score[it], noda_score[it])
    rmse_anl_list.append(rmse_anl)
    rmse_sim_list.append(rmse_sim)

  for it_obs in range(1, obs_timestep):
    obs_true_score = []
    for ix_obs in range(0, nx, obs_xintv):
      obs_true_score.append(true_score[it_obs*obs_tintv][ix_obs])
    rmse_obs = lz.accuracy_rmse_func(obs_true_score, obs_score[it_obs-1])
    rmse_obs_list.append(rmse_obs)

  lz.rmse_draw(timestep, rmse_anl_list, rmse_sim_list, rmse_obs_list, obs_tintv)

  #---------------------------------------------------------- 
  # +++ err cov func.
  #---------------------------------------------------------- 
  """
  for i_num in tqdm(range(0, int(timestep/obs_tintv))):
    matrix_anlta = err.err_anlta[i_num].reshape(nx, nx)
    err.error_heatmap(matrix_anlta, i_num)
  """
  
