#-*- coding: utf-8 -*-
"""
Created on 2020.4.25
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './lib'))

import cal_statics
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D

class lorenz63_score:
  def __init__(self, path:str):
    """Lorenz(1963)モデル結果の描画

    Args:
      path(str): 土台になるデータリストのパス。(lorenz63の並列描画を行いたいデータリスト)
    
    Note:
      self.df(pd.core.frame.DataFrame): 土台データフレーム
      self.time_list(list): 時間リスト
      self.true_list(list): 真値のリスト(OSSEの元)
      self.sim_list(list) : シミュレーション値のリスト
      self.da_list(list)  : データ同化サイクル値のリスト
      
      self.obs(pd.core.frame.DataFrame): 真値リストに乱数を与えて作成した観測値データリスト
      self.obs_list(list): 
      self.obs_time_list(list): 観測値様の時間リスト
      self.obs_true_list(list): 観測値様の時間に相当するtrue_listの値
    """ 
    print('Call ..... Constract of lorenz63_score' + path)
    
    # +++ common list
    self.df  = pd.read_csv(path, sep=',', comment='#', skiprows=2)

    self.time_list = self._df2list(self.df, ' timestep')[0]
    self.true_list = self._df2list(self.df, ' x_true', ' y_true', ' z_true')
    undef = -999.e0
    
    # obs. list
    obs = self.df.loc[ :, [' timestep', ' x_true', ' y_true', ' z_true', ' x_obs', ' y_obs', ' z_obs']]
    obs = obs[obs[' x_obs'] != undef]
    self.obs_time_list = self._df2list(obs, ' timestep')[0]
    self.obs_true_list = self._df2list(obs, ' x_true', ' y_true', ' z_true')
    self.obs_list      = self._df2list(obs, ' x_obs', ' y_obs', ' z_obs')
    
    # sim. list
    self.sim_list  = self._df2list(self.df, ' x_sim', ' y_sim', ' z_sim')
    self.da_list   = self._df2list(self.df, ' x_anl', ' y_anl', ' z_anl')
    
    """
    list.shape
    list[0] = x score, list[1] = y score, list[2] = z score, 
    """
  
  def _df2list(self, Dataframe:pd.core.frame.DataFrame, *index_wrd:tuple) -> list:
    """pandasフレームをリスト化

    Args:
      Dataframe(pd.core.frame.DataFrame): リスト化したいデータフレーム
      index_wrd(tuple): 抽出したいコラム名

    Returns:
      index_list(list): pandasフレームを指定されたコラム名で抽出したリスト
    """
    index_list = []
    for i_index in index_wrd:
      index_list.append(Dataframe[i_index].tolist())
    return index_list



  def lorenz_3ddraw(self, timestep:int) -> int:
    """土台データリストの3D描画

    Args:
        timestep(int): 描画したいタイムステップ

    Returns:
        timestep(int): 忘れた。並列化のtqdmで必要か?
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_xlabel("X Axis"); ax.set_ylabel("Y Axis");  ax.set_zlabel("Z Axis")
    ax.plot(self.true_list[0], self.true_list[1], self.true_list[2], lw=0.5)

    ax.scatter(self.true_list[0][timestep], self.true_list[1][timestep], self.true_list[2][timestep], marker='o', color='b',label='True', alpha=0.5)
    ax.scatter(self.sim_list[0][timestep], self.sim_list[1][timestep], self.sim_list[2][timestep], marker='o', color='r', label='Sim.', alpha=0.5)
    ax.scatter(self.da_list[0][timestep], self.da_list[1][timestep], self.da_list[2][timestep], marker='o', color='g', label='Data assim.', alpha=0.5)

    time = (timestep+1)*0.01
    ax.set_title('Lorenz(1963) - {:.2f} sec.'.format(time), loc='left')
    plt.legend()
    plt.savefig('./figure/Lorenz_xyz_{:0>5}step.png'.format(timestep))
    plt.close('all')

    return timestep

  def lorenz_rmse_draw(self, rmse_sim:list, rmse_obs:list, rmse_kf:list,
                             #rmse_da_1:list, rmse_da_2:list, rmse_da_3:list, rmse_da_4:list, rmse_da_5:list
                             ) -> None:
    """RMSEの描画

    Args:
        rmse_obs(list): 土台データの観測値のRMSE
        rmse_da(list): 土台データのデータ同化のRMSE
        rmse_sim(list): 土台データリストのシミュレーション値のRMSE
        +
        *適宜、加える。(可変引数を使うことはできるか??)
        rmse_da_?(list): 加えたい違う種類のデータ同化のRMSE
    """
    fig = plt.figure()
    ax1 = fig.subplots()
    
    sns.set_style('whitegrid')
    ax1.plot(self.time_list, rmse_sim, ls="--", color='b', label='SIM')
    ax1.plot(self.time_list, rmse_kf, ls="--", color='r', label='KF')
    
    #ax1.plot(self.time_list, rmse_da_1, ls=":", alpha=0.3, label='EnKF 3mem')
    #ax1.plot(self.time_list, rmse_da_2, ls="--", label='EnKF 10mem')
    #ax1.plot(self.time_list, rmse_da_3, ls="--", label='EnKF 100mem')
    #ax1.plot(self.time_list, rmse_da_4, ls="--", label='EnKF 1000mem')
    #ax1.plot(self.time_list, rmse_da_5, ls="--", label='EnKF 5000mem')

    ax1.scatter(self.obs_time_list, rmse_obs, marker='*', color='y', s=35, alpha=0.5, edgecolor='k', label='Obs.')

    ax1.set_xlabel('sec.')
    ax1.set_ylabel('RMSE')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 10)
    ax1.set_title('Lorenz(1963) RMSE', loc='left')

    plt.grid()
    plt.legend()
    plt.show()

class lorenz63_errcov:
  """
  一つのデータリストを誤差共分散行列の扱いを行うクラスメソッド
  """
  def __init__(self, path:str, *, matrix_size:int=3):
    self.err_data = self._read_errcov(path)
    
  def _read_errcov(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")

  def trace_errcov(self, matrix) -> np.ndarray:
    return np.trace(matrix)

  def error_heatmap(self, err_data:np.ndarray, timestep:int) -> None:
    ax = plt.subplots()
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.set(style='white')
    sns.heatmap(
      data=err_data, cmap=cmap, annot=True, fmt='.6f',
      vmax=1.0, vmin=-1.0, center=0,
      square=True, linewidths=.5, 
      cbar_kws={"shrink": .5, "extend": 'both'},
      xticklabels=['x','y','z'], yticklabels=['x','y','z']
    )

    time = timestep*0.01
    ax.set_title('Back ground Error Cov. - {:.1f} sec.'.format(time), loc='center')
    plt.savefig('./figure/error_heatmap_{:.1f}sec.png'.format(time))
    plt.close('all')

""" private package """

def accuracy_rmse_func(true_list:list, asses_list:list, num_elem:int) -> list:
  """RMSEリストの作成
  Args:
      true_list(list): 今回の実験の真値
      asses_list(list): 真値との差を図りたいリスト
      num_elem(int): 真値との差を図りたい変数の数
  Returns:
      rmse_list(list): タイムステップ別に分けたrmse
  Note:
    シミュレーション値, データ同化値と観測値では, RMSEのタイムリストや変数の数が違うので注意。
  """
  rmse_list  = []
  time_range = len(true_list[0])
  
  for i_num in range(time_range):
    rmse_true_list , rmse_asses_list= [], []
    for i_elem in range(num_elem):
      rmse_true_list.append(true_list[i_elem][i_num])
      rmse_asses_list.append(asses_list[i_elem][i_num])
    rmse = cal_statics.rmse(rmse_true_list, rmse_asses_list)
    rmse_list.append(rmse)
  
  return rmse_list

def _enkf_pathset(mem:int, method_num:int, *, outdir:str='./output/lorenz63/') -> str:
  """EnKFのPATH設定

  Args:
      mem(int): EnKFのメンバー数
      method_num(int): 0はPO法、1はSRF法

  Keyword Arguments:
      outdir(str): PATHの大元設定。変更なければ引数に加えなくて良い。 (default: {'./output/lorenz63'})

  Returns:
      enkf_data_path(str):  enkfのデータPATH
      enkf_errcov_path(str):  enkfの誤差共分散データPATH
  """
  enkf_method = ['PO', 'SRF']
  enkf_data_path = outdir+'EnKF_'+str(mem)+'m_'+enkf_method[method_num]+'.csv'
  enkf_errcov_path = outdir+'errcov_EnKF_'+str(mem)+'m_'+enkf_method[method_num]+'.csv'
  return enkf_data_path, enkf_errcov_path

if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  #---------------------------------------------------------- 
  matrix_size  = 3
  obs_interval = 1

  outdir = './output/lorenz63/'
  path_kf, path_kf_errcov = outdir+'KF.csv', outdir+'errcov_KF.csv'
  
  """Ensemble data set"""
  # for compare inflation
  path_enkf_2m_PO, path_enkf_2m_errcov_PO = _enkf_pathset(2,0)
  path_enkf_2m_SRF, path_enkf_2m_errcov_SRF = _enkf_pathset(2,1)
  # for compare member
  path_enkf_5m_SRF, path_enkf_5m_errcov_SRF = _enkf_pathset(5,1)
  path_enkf_15m_SRF, path_enkf_15m_errcov_SRF = _enkf_pathset(15,1)
  path_enkf_50m_SRF, path_enkf_50m_errcov_SRF = _enkf_pathset(50,1)

  #---------------------------------------------------------- 
  # +++ class set
  # >>  Args (1) lorenz63 result. score
  # >>       (2) Pf (rediction err covariance matrix)
  #---------------------------------------------------------- 
  score_kf, errcov_kf = lorenz63_score(path_kf), lorenz63_errcov(path_kf_errcov)
  score_enkf_2m_PO, errcov_enkf_2m_PO = lorenz63_score(path_enkf_2m_PO), lorenz63_errcov(path_enkf_2m_errcov_PO)
  score_enkf_2m_SRF, errcov_enkf_2m_SRF = lorenz63_score(path_enkf_2m_SRF), lorenz63_errcov(path_enkf_2m_errcov_SRF)
  
  score_enkf_5m_SRF = lorenz63_score(path_enkf_5m_SRF)
  score_enkf_15m_SRF = lorenz63_score(path_enkf_15m_SRF)
  score_enkf_50m_SRF = lorenz63_score(path_enkf_50m_SRF)

  #---------------------------------------------------------- 
  # +++ draw func.
  # >>  (1) RMSE draw
  # >>  (2) 3D draw
  #---------------------------------------------------------- 

  """ (1) RMSE draw """
  num_sim_elem, num_obs_elem = 3, 3
  rmse_sim = accuracy_rmse_func(score_kf.true_list, score_kf.sim_list, num_sim_elem)
  rmse_kf  = accuracy_rmse_func(score_kf.true_list, score_kf.da_list, num_sim_elem)
  rmse_obs = accuracy_rmse_func(score_kf.obs_true_list, score_kf.obs_list, num_obs_elem)

  score_kf.lorenz_rmse_draw(rmse_sim, rmse_obs, rmse_kf,)

  sys.exit()

  """ (2) 3D draw @using Pool(parallel method) """
  laststep = 2500
  viewstep = int(laststep/3)+1
  time_before = time.time()
  values = [ i_num for i_num in range(0, len(score_kf.time_list[0:laststep]), 3)]

  with Pool(4) as pool:
    sleeps = list(tqdm(pool.imap(score_kf.lorenz_3ddraw, values), total=viewstep))

  time_after = time.time()
  elapsed_time = time_after - time_before
  print(f'time: {elapsed_time} sec.')