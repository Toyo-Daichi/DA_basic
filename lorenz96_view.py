#-*- coding: utf-8 -*-
"""
Created on 2020.5.7
@author: Toyo_anlichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './lib'))

import cal_statics
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

class lorenz96_score:
  def __init__(self):
    print('..... Constract of lorenz96_score class ')

  def lorenz_rmse_draw(self, timestep:int, rmse_anl:list, rmse_sim:list, rmse_obs:list, obs_tintv:int) -> None:
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

  def lorenz_hov_draw(self, nx:int, timestep:int, lorenz_anlta:np.ndarray) -> None:
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
    self.err_anlta = _csv2list(path)

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

""" private package """

def _csv2list(path:str) -> np.ndarray:
  return np.genfromtxt(path, delimiter=",")

def _accuracy_rmse_func(true_list:list, asses_list:list) -> list:
  rmse = cal_statics.rmse(true_list, asses_list)
  return rmse

if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  #---------------------------------------------------------- 
  # (1) dimension
  obs_xintv = 1
  obs_tintv = 1
  day_tintv = 4
  init_step = 1
  nx = 40
  ny = int(nx/obs_xintv)

  stepday = 40
  timestep = stepday*day_tintv

  timeshape = init_step + timestep
  obs_timeshape =  int(timestep/obs_tintv)

  # (2) PATH setting
  # >> The three file names below do not change in any way.
  outdir = './output/lorenz96/'
  true_path = outdir + 'normal_true_score_' + str(nx) + 'ndim.csv'
  sim_path = outdir + 'normal_sim_score_' + str(nx) + 'ndim.csv'
  obs_path = outdir + 'normal_obs_score_' + str(nx) + 'ndim.csv'

  """ KF data set """
  path_kf = outdir + 'normal_KF_anl_score_' + str(nx) + 'ndim.csv'
  path_kf_errcov = outdir + 'normal_KF_anlinc_' + str(nx) + 'ndim.csv'
  path_kf_anlinc = outdir + 'normal_KF_anlinc_' + str(nx) + 'ndim.csv'
  
  """Ensemble KF data set (basic is EnSRF.) """
  mems = 50
  path_enkf = outdir + 'normal_EnKF' + str(mems) + 'm_anl_score_' + str(nx) + 'ndim.csv'
  path_kf_errcov = outdir + 'normal_EnKF' + str(mems) + 'm_anlinc_' + str(nx) + 'ndim.csv'
  path_kf_anlinc = outdir + 'normal_EnKF' + str(mems) + 'm_errcov_' + str(nx) + 'ndim.csv'

  #---------------------------------------------------------- 
  # +++ reading basic score.
  #---------------------------------------------------------- 
  true_score = _csv2list(true_path).reshape(timeshape, nx)
  sim_score = _csv2list(sim_path).reshape(timeshape, nx)
  obs_score = _csv2list(obs_path).reshape(obs_timeshape, ny)

  #---------------------------------------------------------- 
  # +++ class set
  # (1) lorenz96_score
  #---------------------------------------------------------- 
  lorenz96 = lorenz96_score()
  
  kf_anl_score = _csv2list(path_kf).reshape(timeshape, nx)
  enkf_anl_score = _csv2list(path_enkf).reshape(timeshape, nx)

  #---------------------------------------------------------- 
  # +++ RMSE func.
  #---------------------------------------------------------- 
  