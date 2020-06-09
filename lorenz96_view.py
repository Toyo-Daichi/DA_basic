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
import statics_tool
from tqdm import tqdm
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

class lorenz96_score:
  def __init__(self):
    pass

  def trajectory_draw(
    self, true_score:np.ndarray, sim_score:np.ndarray, anl_score:np.ndarray, obs_score:np.ndarray,
    *, time_length:int=70, obs_tintv:int=1
    ) -> None:
    """一つの要素の時間変化を示す図の作成
    Args:
        true_score (np.ndarray): 真値
        sim_score (np.ndarray): シミュレーション値
        anl_score (np.ndarray): データ同化解析値(***今後、必要になればこのshapeを増やす)
        obs_score (np.ndarray): 真値にランダムノイズを加えた観測値
    Note:
        Shape of true_score, sim_score, anl_score is [timescale, nx].
        Shape of obs_score is [obs_timescale, ny]. OBS dont have 0step score.
    """
    fig, ax = plt.subplots(2, 3, figsize=(15,10))
    x_coord = np.arange(time_length)
    obs_x_coord = np.arange(1,time_length, obs_tintv)

    self._trajectory_basic(0,0,0,time_length,ax,x_coord,obs_x_coord,true_score,sim_score,anl_score,obs_score)
    self._trajectory_basic(0,1,1,time_length,ax,x_coord,obs_x_coord,true_score,sim_score,anl_score,obs_score)
    self._trajectory_basic(0,2,2,time_length,ax,x_coord,obs_x_coord,true_score,sim_score,anl_score,obs_score)
    self._trajectory_basic(1,0,3,time_length,ax,x_coord,obs_x_coord,true_score,sim_score,anl_score,obs_score)
    self._trajectory_basic(1,1,4,time_length,ax,x_coord,obs_x_coord,true_score,sim_score,anl_score,obs_score)
    self._trajectory_basic(1,2,5,time_length,ax,x_coord,obs_x_coord,true_score,sim_score,anl_score,obs_score)

    plt.savefig('./figure/trajectory.png')

  def _trajectory_basic(
    self, fig_line:int, fig_column:int, xdim:int, time_length:int,
    ax, x_coord:np.ndarray, obs_x_coord:np.ndarray,
    true_score:np.ndarray, sim_score:np.ndarray, anl_score:np.ndarray, obs_score:np.ndarray
    ) -> None:
    """trajectory_draw図の基本情報
    Args:
        fig_line, fig_colum (int): 図の配置場所
        xdim (int): 書きたい変数番号
        time_length (int): 時間レンジの指定
    """
    _draw=ax[fig_line, fig_column]
    _draw.set_xlabel('TIME STEP')
    _draw.set_ylabel('X ' + str(xdim) + ' dim. COORDINATE')
    _draw.set_xlim(0, time_length)
    _draw.set_title('LORENZ(1996) X ' + str(xdim) + ' dim. COORDINATE', loc='left')
    _draw.plot(x_coord, true_score[0:time_length, xdim], lw=2.0, ls="--", color='b', label='TRUE')
    _draw.plot(x_coord, anl_score[0:time_length, xdim], ls="-", color='r', label='DATA ASSIM')
    _draw.scatter(obs_x_coord, obs_score[0:time_length-1, xdim], marker='*',  s=35, alpha=0.5, edgecolor='k', label='OBS')
    _draw.grid()
    _draw.legend(loc='upper right')

  def hovmeller_draw(
    self, true_score:np.ndarray, sim_score:np.ndarray, anl_score:np.ndarray, anl_inc_score:np.ndarray,
    *, nx=40, timeshape=161) -> None:
    """それぞれのデータセットでのホフメラ図
    Args:
        true_score (np.ndarray): 真値
        sim_score (np.ndarray): シミュレーション値
        anl_score (np.ndarray): データ同化解析値
        anl_inc_score (np.ndarray): 解析インクリメント値
    """
    fig, ax = plt.subplots(1,3, figsize=(20,5))
    x_coord = np.arange(nx)
    y_coord, y_coord_inc = np.arange(timeshape), np.arange(timeshape-1)
   
    _, _ = self._hovmeller_basic(0,50,true_score,'TRUE',ax,x_coord,y_coord)
    _, _ = self._hovmeller_basic(1,50,sim_score,'SIMULATION',ax,x_coord,y_coord)
    _ax, _cbar = self._hovmeller_basic(2,50,anl_score,'DATA ASSIM',ax,x_coord,y_coord)
    fig.colorbar(_cbar, ax=_ax)
    plt.savefig('./figure/houvmeller.png')

    fig, ax = plt.subplots(1,2, figsize=(11,5))
    levels = np.arange(-1.0, 1.0, 0.01)
    _ax1, _cbar1 = self._hovmeller_basic(0,70,true_score-sim_score,'TRUE - SIMULATION',ax,x_coord,y_coord, cmap=plt.get_cmap('bwr'))
    _ax2, _cbar2 = self._hovmeller_basic(1,70,true_score-anl_score,'TRUE - DATA ASSIM',ax,x_coord,y_coord, cmap=plt.get_cmap('bwr'), levels=levels)
    fig.colorbar(_cbar1, ax=_ax1)
    fig.colorbar(_cbar2, ax=_ax2)
    plt.savefig('./figure/houvmeller_diff.png')

    fig, ax = plt.subplots(1,1, figsize=(6,5))
    _ax, _cbar = self._hovmeller_inc(50,anl_inc_score,'DATA ASSIM INCREMENT',ax,x_coord,y_coord_inc)
    fig.colorbar(_cbar, ax=_ax)
    plt.savefig('./figure/houvmeller_inc.png')

  def _hovmeller_basic(
    self, fig_line:int, time_length:int, data:np.ndarray, data_form:str, 
    ax, x_coord:np.ndarray, y_coord:np.ndarray,
    *, cmap=plt.get_cmap('PiYG'), levels=np.arange(-8.0,8.0,0.1)
    ) -> None:
    """hovmeller_draw図の基本情報
    Args:
        fig_line (int): 図の配置場所
        data (np.ndarray), data_form(str): 書きたいデータセット, データセットの命名
        time_length (int): 時間レンジの指定
    """
    _draw = ax[fig_line]
    _cbar = _draw.contourf(x_coord, y_coord, data, levels, cmap=cmap, extend='both')
    _draw.set_xlabel('X COORDINATE')
    _draw.set_ylabel('TIME STEP')
    _draw.set_ylim(time_length,0)
    _draw.set_xlim(nx-1,0)
    _draw.set_title('LORENZ(1996) ' + data_form, loc='left')

    return _draw, _cbar

  def _hovmeller_inc(
    self, time_length:int, data:np.ndarray, data_form:str, 
    ax, x_coord:np.ndarray, y_coord:np.ndarray) -> None:
    """hovmeller_draw図の解析インクリメントの時間変化
    Args:
        fig_line (int): 図の配置場所
        data (np.ndarray), data_form(str): 書きたいデータセット, データセットの命名
        time_length (int): 時間レンジの指定
    Note:
        Shape of anl_inc_score is [obs_timescale, nx]. ANL_INC dont have 0step score.
    """
    levels, cmap = np.arange(-8.0, 8.0, 0.05), plt.get_cmap('bwr')
    _draw = ax
    _cbar = _draw.contourf(x_coord, y_coord, data, levels, cmap=cmap, extend='both')
    _draw.set_xlabel('X COORDINATE')
    _draw.set_ylabel('TIME STEP')
    _draw.set_ylim(time_length,0)
    _draw.set_xlim(nx-1,0)
    _draw.set_title('LORENZ(1996) ' + data_form, loc='left')

    return _draw, _cbar
    
  def rmse_draw(
    self, rmse_sim_list:list, rmse_obs_list:list, rmse_anl_list:list,
    # add other kind rmse_anl_list
    rmse_v1_enkf_anl_list,
    rmse_v2_enkf_anl_list,
    rmse_v3_enkf_anl_list,
    rmse_v4_enkf_anl_list,
    rmse_v5_enkf_anl_list,
    *, time_length=50, timeshape=161, obs_timeshape=160
    ) -> None:
    """LORENZ(1996)での双子実験で作成した真値とのRMSE
    Args:
        rmse_sim_list (list): 真値からシミュレーション値を引いて作成した値のリスト
        rmse_obs_list (list): 真値から観測値を引いて作成した値のリスト
        rmse_anl_list (list): 真値から解析値を引いて作成した値のリスト

        time_length (int, optional): 描きたいタイムのレンジの幅。 Defaults to 50.
        timeshape (int, optional): データが所持しているタイムステップ幅。 Defaults to 161.
        obs_timeshape (int, optional): 観測データが所持しているタイムステップ幅。 Defaults to 160.
    """
    fig, ax = plt.subplots(figsize=(6, 5))
    time_list, obs_time_list = np.arange(timeshape), np.arange(obs_timeshape)

    _draw = ax
    _draw.plot(time_list[0:time_length], rmse_sim_list[0:time_length], ls=":", color='b', label='SIMULATION')
    _draw.plot(time_list[0:time_length], rmse_anl_list[0:time_length], ls="-", color='r', label='DATA ASSIM EKF')

    # for added other kind rmse_anl_list
    _draw.plot(time_list[0:time_length], rmse_v1_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=20m')
    #_draw.plot(time_list[0:time_length], rmse_v1_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=20m@loc')
    _draw.plot(time_list[0:time_length], rmse_v2_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=40m')
    #_draw.plot(time_list[0:time_length], rmse_v2_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=40m@loc')
    _draw.plot(time_list[0:time_length], rmse_v3_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=100m')
    #_draw.plot(time_list[0:time_length], rmse_v3_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=100m@loc')
    _draw.plot(time_list[0:time_length], rmse_v4_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=500m')
    _draw.plot(time_list[0:time_length], rmse_v5_enkf_anl_list[0:time_length], ls="--", label='DATA ASSIM EnKF=1000m')

    _draw.scatter(obs_time_list[0:time_length], rmse_obs_list[0:time_length], marker='*', color='y', s=35, alpha=0.5, edgecolor='k', label='OBS')
    _draw.set_xlabel('TIMESTEP')
    _draw.set_ylabel('RMSE')
    _draw.set_title('LORENZ(1996) RMSE', loc='left')
    _draw.set_ylim(0,7)
    _draw.grid()
    _draw.legend(loc='upper right')
    plt.savefig('./figure/lorenz96.png')

  def making_rmse_snap(
    self, true_score, sim_score, anl_score, obs_score, timeshape, obs_timeshape, obs_xintv
    ) -> list:
    """1ステップにあたる真値とのRMSEリストの作成
    Args:
        true_score (np.ndarray): 真値
        sim_score (np.ndarray): シミュレーション値
        anl_score (np.ndarray): データ同化解析値
        obs_score (np.ndarray): 真値にランダムノイズを加えた観測値

        timeshape (int): データが持つタイムステップ数
        obs_timeshape (int): 観測データが持つタイムステップ数(0ステップは持っていない)
        obs_xintv (int) : 観測データがどの真値にノイズを加えたかの判別に必要

    Returns:
        list: 各リストの時間ごとのRMSEを求めたリスト
    """
    
    rmse_anl_list, rmse_sim_list, rmse_obs_list = [], [], []

    for _it in range(timeshape):
      rmse_anl = _accuracy_rmse_func(true_score[_it], anl_score[_it])
      rmse_sim = _accuracy_rmse_func(true_score[_it], sim_score[_it])
      rmse_anl_list.append(rmse_anl)
      rmse_sim_list.append(rmse_sim)

    for _it in range(obs_timeshape):
      obs_true_score = []
      for _obs in range(0, nx, obs_xintv):
        obs_true_score.append(true_score[_it*obs_tintv+obs_tintv][_obs])
      rmse_obs = _accuracy_rmse_func(obs_true_score, obs_score[_it])
      rmse_obs_list.append(rmse_obs)

    return rmse_anl_list, rmse_sim_list, rmse_obs_list

class lorenz96_errcov:
  def __init__(self):
    pass

  def errcov_draw(self, errcov_mtx:np.ndarray, *, nx=40, state='Pa', time=1) -> None:
    """誤差共分散行列の描画
    Args:
        errcov_mtx (np.ndarray): 時間が指定してある誤差共分散-np.ndarray
        nx (int) : 変数の数
    """
    fig, ax = plt.subplots()
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.set(style='white')
    sns.heatmap(
      data=errcov_mtx, cmap=cmap, 
      vmax=0.1, vmin=-0.1, center=0,
      square=True, linewidths=.5, 
      cbar_kws={"shrink": .5, "extend": 'both'},
      xticklabels=2, yticklabels=2
    )
    ax.set_title('LORENZ(1996) ERROR COVARIANCE MATRIX ({}), {:d} step.'.format(state, time), fontsize=8, loc='left')
    plt.yticks(rotation=0)
    plt.savefig('./figure/errcov_matrix_{:03d}step.png'.format(time))
    plt.close('all')

  def cross_corr_draw(self, cross_corr, cross_cov, y_gauss, *, nx=40, time=1) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    _cor = ax
    _cov = ax.twinx()
    x_coord = np.arange(nx)

    _cor.scatter(x_coord, cross_corr, marker='*', color='y')
    lns1 = _cor.plot(x_coord, cross_corr, ls='-', color='r', label='CORRELATION')
    lns2 = _cov.plot(x_coord, cross_cov, ls='-', color='k', lw=2, alpha= 0.5, label='COVARIANCE')
    lns3 = _cor.plot(x_coord, y_gauss, ls=':', color='b', alpha=0.5, label='GAUSSIAN FUNCTION')
    
    _cor.set_xlabel('X COORDINATE')
    _cor.set_ylabel('CORRELATION')
    _cov.set_ylabel('COVARIANCE')
    
    _cor.set_ylim(-1.0,1.0)
    _cov.set_ylim(-0.3,0.3)

    _cor.set_title('LORENZ(1996) CROSS CORRELATION {:d} step.'.format(time), loc='left')
    # added these three lines
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    _cor.legend(lns, labs, loc=0)
    plt.savefig('./figure/errcov_cross_{:03d}step.png'.format(time))
    plt.close('all')

  def making_cross_corr(self, anl_errcov:np.ndarray, obs_timeshape:int, *, target_grd=19) -> list:
    """相関係数の断面図的なlistの作成??
    Args:
        anl_errcov (np.ndarray): 誤差共分散データ(時系列含む)
        obs_timeshape (int): obs_timeshape (int): 観測データが持つタイムステップ数(0ステップは持っていない)
        target_grd (int, optional): クロス断面図を作るにあたる相関をとるターゲットのグリッド。Defaults to 20.
    Returns:
        list: target_grdとの各地点での相関リスト(時系列含む)
    """
    cross_corr = [[]*i for i in range(obs_timeshape)]
    gaussian_list = [[]*i for i in range(obs_timeshape)]

    for _it in range(obs_timeshape):
      diag = np.diag(anl_errcov[_it])
      errcov_line = anl_errcov[_it, target_grd]
    
      for _ix in range(nx):
        corr = errcov_line[_ix]/(math.sqrt(diag[_ix])*math.sqrt(diag[target_grd]))
        cross_corr[_it].append(corr)
    
      gaussian_list[_it] = self._compare_gaussian(cross_corr[_it])
    return cross_corr, gaussian_list

  def _compare_gaussian(self, cross_corr:list, *, nx=40, target_grd=19) -> list:
    gaussian_list = []
    std_list = np.arange(0.05, 2.0, 0.05)
    x_coord = np.arange(nx)
    for _ in std_list:
      y_gauss = _gaussian_func(x_coord, amp=1.0, ave=target_grd, std=_)
      gaussian_list.append(y_gauss)
    
    rmse_list = []
    for _ in gaussian_list:
      rmse_list.append(_accuracy_rmse_func(cross_corr, _))
    
    least_case = np.argmin(rmse_list)
    #least_rmse = min(rmse_list)
    return gaussian_list[least_case]


""" private package """

def _gaussian_func(x, *, amp:float=1.0, ave:float=0, std:float=1.0):
    """ガウス関数
    Args:
        amp(float): 振幅
        ave(float): 平均
        std(float): 標準偏差
    Returns:
        gaussian function(float)
    """
    return amp*np.exp(-1*((x - ave)/2*std)**2)

def _csv2list(path:str) -> np.ndarray:
  return np.genfromtxt(path, delimiter=",")

def _accuracy_rmse_func(true_list:list, asses_list:list) -> list:
  rmse = statics_tool.rmse(true_list, asses_list)
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

  #Notes: OBS dont have 0step score.
  timeshape = init_step + timestep
  obs_timeshape = int(timestep/obs_tintv)

  # (2) PATH setting
  # >> The three file names below do not change in any way.
  outdir = './output/lorenz96/'
  true_path = outdir + 'normal_true_score_' + str(nx) + 'ndim.csv'
  sim_path = outdir + 'normal_sim_score_' + str(nx) + 'ndim.csv'
  obs_path = outdir + 'normal_obs_score_' + str(nx) + 'ndim.csv'

  """ KF data set """
  path_kf = outdir + 'normal_KF_anl_score_' + str(nx) + 'ndim.csv'
  path_kf_errcov = outdir + 'normal_KF_errcov_' + str(nx) + 'ndim.csv'
  path_kf_anlinc = outdir + 'normal_KF_anlinc_' + str(nx) + 'ndim.csv'
  
  """Ensemble KF data set (basic is EnSRF.) """
  mems, loc = 20, ''
  v1_path_enkf = outdir + 'normal_EnKF' + str(mems) + 'm_anl_score_' + str(nx) + 'ndim' + loc +'.csv'
  v1_path_enkf_errcov = outdir + 'normal_EnKF' + str(mems) + 'm_errcov_' + str(nx) + 'ndim' + loc +'.csv'
  v1_path_enkf_anlinc = outdir + 'normal_EnKF' + str(mems) + 'm_anlinc_' + str(nx) + 'ndim' + loc +'.csv'
  mems, loc = 40, ''
  v2_path_enkf = outdir + 'normal_EnKF' + str(mems) + 'm_anl_score_' + str(nx) + 'ndim' + loc +'.csv'
  v2_path_enkf_errcov = outdir + 'normal_EnKF' + str(mems) + 'm_errcov_' + str(nx) + 'ndim' + loc +'.csv'
  v2_path_enkf_anlinc = outdir + 'normal_EnKF' + str(mems) + 'm_anlinc_' + str(nx) + 'ndim' + loc +'.csv'
  mems, loc = 100, ''
  v3_path_enkf = outdir + 'normal_EnKF' + str(mems) + 'm_anl_score_' + str(nx) + 'ndim' + loc +'.csv'
  v3_path_enkf_errcov = outdir + 'normal_EnKF' + str(mems) + 'm_errcov_' + str(nx) + 'ndim' + loc +'.csv'
  v3_path_enkf_anlinc = outdir + 'normal_EnKF' + str(mems) + 'm_anlinc_' + str(nx) + 'ndim' + loc +'.csv'
  mems, loc = 500, ''
  v4_path_enkf = outdir + 'normal_EnKF' + str(mems) + 'm_anl_score_' + str(nx) + 'ndim' + loc +'.csv'
  v4_path_enkf_errcov = outdir + 'normal_EnKF' + str(mems) + 'm_errcov_' + str(nx) + 'ndim' + loc +'.csv'
  v4_path_enkf_anlinc = outdir + 'normal_EnKF' + str(mems) + 'm_anlinc_' + str(nx) + 'ndim' + loc +'.csv'
  mems, loc = 1000, ''
  v5_path_enkf = outdir + 'normal_EnKF' + str(mems) + 'm_anl_score_' + str(nx) + 'ndim' + loc +'.csv'
  v5_path_enkf_errcov = outdir + 'normal_EnKF' + str(mems) + 'm_errcov_' + str(nx) + 'ndim' + loc +'.csv'
  v5_path_enkf_anlinc = outdir + 'normal_EnKF' + str(mems) + 'm_anlinc_' + str(nx) + 'ndim' + loc +'.csv'

  #---------------------------------------------------------- 
  # +++ reading basic score.
  #---------------------------------------------------------- 
  true_score = _csv2list(true_path).reshape(timeshape, nx)
  sim_score = _csv2list(sim_path).reshape(timeshape, nx)
  obs_score = _csv2list(obs_path).reshape(obs_timeshape, ny)

  #---------------------------------------------------------- 
  # +++ class set
  #  (1) lorenz96_score
  #  >> & anl_data set
  #---------------------------------------------------------- 
  """
  lorenz96_score = lorenz96_score()
  
  kf_anl_score = _csv2list(path_kf).reshape(timeshape, nx)
  kf_anlinc_score = _csv2list(path_kf_anlinc).reshape(obs_timeshape, nx)

  v1_enkf_anl_score = _csv2list(v1_path_enkf).reshape(timeshape, nx)
  v1_enkf_anlinc_score = _csv2list(v1_path_enkf_anlinc).reshape(obs_timeshape, nx)

  v2_enkf_anl_score = _csv2list(v2_path_enkf).reshape(timeshape, nx)
  v2_enkf_anlinc_score = _csv2list(v2_path_enkf_anlinc).reshape(obs_timeshape, nx)

  v3_enkf_anl_score = _csv2list(v3_path_enkf).reshape(timeshape, nx)
  v3_enkf_anlinc_score = _csv2list(v3_path_enkf_anlinc).reshape(obs_timeshape, nx)

  v4_enkf_anl_score = _csv2list(v4_path_enkf).reshape(timeshape, nx)
  v4_enkf_anlinc_score = _csv2list(v4_path_enkf_anlinc).reshape(obs_timeshape, nx)

  v5_enkf_anl_score = _csv2list(v5_path_enkf).reshape(timeshape, nx)
  v5_enkf_anlinc_score = _csv2list(v5_path_enkf_anlinc).reshape(obs_timeshape, nx)
  """

  #---------------------------------------------------------- 
  # +++ Trajectory & Hovmeller 
  #---------------------------------------------------------- 
  #lorenz96_score.trajectory_draw(true_score, sim_score, kf_anl_score, obs_score)
  #lorenz96_score.hovmeller_draw(true_score, sim_score, kf_anl_score, kf_anlinc_score)

  #---------------------------------------------------------- 
  # +++ RMSE 
  #---------------------------------------------------------- 
  """
  rmse_anl_list,rmse_sim_list,rmse_obs_list = \
  lorenz96_score.making_rmse_snap(true_score, sim_score, kf_anl_score, obs_score, timeshape, obs_timeshape, obs_xintv)
  rmse_v1_enkf_anl_list, _, _ = \
  lorenz96_score.making_rmse_snap(true_score, sim_score, v1_enkf_anl_score, obs_score, timeshape, obs_timeshape, obs_xintv)
  rmse_v2_enkf_anl_list, _, _ = \
  lorenz96_score.making_rmse_snap(true_score, sim_score, v2_enkf_anl_score, obs_score, timeshape, obs_timeshape, obs_xintv)
  rmse_v3_enkf_anl_list, _, _ = \
  lorenz96_score.making_rmse_snap(true_score, sim_score, v3_enkf_anl_score, obs_score, timeshape, obs_timeshape, obs_xintv)
  rmse_v4_enkf_anl_list, _, _ = \
  lorenz96_score.making_rmse_snap(true_score, sim_score, v4_enkf_anl_score, obs_score, timeshape, obs_timeshape, obs_xintv)
  rmse_v5_enkf_anl_list, _, _ = \
  lorenz96_score.making_rmse_snap(true_score, sim_score, v5_enkf_anl_score, obs_score, timeshape, obs_timeshape, obs_xintv)
  """

  #lorenz96_score.rmse_draw(rmse_sim_list, rmse_obs_list, rmse_anl_list,
  #rmse_v1_enkf_anl_list, rmse_v2_enkf_anl_list, rmse_v3_enkf_anl_list, rmse_v4_enkf_anl_list, rmse_v5_enkf_anl_list,
  #time_length=100
  #)

  #---------------------------------------------------------- 
  # +++ class set
  #  (2) lorenz96_errcov
  #  >> & anl_data set
  #---------------------------------------------------------- 
  lorenz96_errcov = lorenz96_errcov()
  kf_anl_errcov = _csv2list(path_kf_errcov).reshape(obs_timeshape,nx,nx)
  v1_enkf_anl_errcov = _csv2list(v1_path_enkf_errcov).reshape(obs_timeshape,nx,nx)
  v2_enkf_anl_errcov = _csv2list(v2_path_enkf_errcov).reshape(obs_timeshape,nx,nx)
  v3_enkf_anl_errcov = _csv2list(v3_path_enkf_errcov).reshape(obs_timeshape,nx,nx)
  v4_enkf_anl_errcov = _csv2list(v4_path_enkf_errcov).reshape(obs_timeshape,nx,nx)
  v5_enkf_anl_errcov = _csv2list(v5_path_enkf_errcov).reshape(obs_timeshape,nx,nx)

  filter_stop_step = 100
  for _it in tqdm(range(100)):
    #lorenz96_errcov.errcov_draw(kf_anl_errcov[_it], time=_it)
    lorenz96_errcov.errcov_draw(v1_enkf_anl_errcov[_it], state='Pf', time=_it)
    #lorenz96_errcov.errcov_draw(v2_enkf_anl_errcov[_it], state='Pf', time=_it)
    #lorenz96_errcov.errcov_draw(v3_enkf_anl_errcov[_it], state='Pf', time=_it)
    #lorenz96_errcov.errcov_draw(v4_enkf_anl_errcov[_it], state='Pf', time=_it)
    #lorenz96_errcov.errcov_draw(v5_enkf_anl_errcov[_it], state='Pf', time=_it)

  #---------------------------------------------------------- 
  # +++ Cross correlation list
  #---------------------------------------------------------- 
  #cross_corr_list, gaussian_list = lorenz96_errcov.making_cross_corr(kf_anl_errcov, obs_timeshape)
  #for _it in tqdm(range(obs_timeshape)):
  #  cross_cov_list = kf_anl_errcov[_it, 19]
  #  lorenz96_errcov.cross_corr_draw(cross_corr_list[_it], cross_cov_list, gaussian_list[_it], time=_it+1)

  print('Normal END')






  