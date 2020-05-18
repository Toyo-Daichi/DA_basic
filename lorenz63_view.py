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
    print('Call ..... Constract of lorenz63_score')
    
    # +++ common list
    self.df  = pd.read_csv(path, sep=',')

    self.time_list = self.df2_list(self.df, ' timestep')[0]
    self.true_list = self.df2_list(self.df, ' x_true', ' y_true', ' z_true')
    undef = -999.e0
    
    # obs. list
    obs = self.df.loc[ :, [' timestep', ' x_true', ' y_true', ' x_obs', ' y_obs']]
    obs = obs[obs[' x_obs'] != undef]
    self.obs_time_list  = self.df2_list(obs, ' timestep')[0]
    self.obs_true_list = self.df2_list(obs, ' x_true', ' y_true')
    self.obs_list       = self.df2_list(obs, ' x_obs', ' y_obs')
    
    # sim. list
    self.sim_list  = self.df2_list(self.df, ' x_sim', ' y_sim', ' z_sim')
    self.da_list   = self.df2_list(self.df, ' x_da', ' y_da', ' z_da')
    
    """
    list.shape
    list[0] = x score, list[1] = y score, list[2] = z score, 
    """
  
  def df2_list(self, Dataframe:pd.core.frame.DataFrame, *index_wrd:tuple) -> list:
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

  def accuracy_rmse_func(self, true_list:list, asses_list:list, num_elem:int) -> list:
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

  def lorenz_rmse_draw(self, rmse_sim:list, rmse_obs:list, rmse_da:list,
                             #rmse_da_1:list, rmse_da_2:list, rmse_da_3:list, rmse_da_4:list, rmse_da_5:list
                             ) -> None:
    """RMSEの描画

    Arguments:
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
    ax1.plot(self.time_list, rmse_sim, ls="--", color='b', label='No DA')
    ax1.plot(self.time_list, rmse_da, ls="--", color='r', label='KF')
    
    #ax1.plot(self.time_list, rmse_da_1, ls=":", alpha=0.3, label='EnKF 3mem')
    #ax1.plot(self.time_list, rmse_da_2, ls="--", label='EnKF 10mem')
    #ax1.plot(self.time_list, rmse_da_3, ls="--", label='EnKF 100mem')
    #ax1.plot(self.time_list, rmse_da_4, ls="--", label='EnKF 1000mem')
    #ax1.plot(self.time_list, rmse_da_5, ls="--", label='EnKF 5000mem')

    ax1.scatter(self.obs_time_list, rmse_obs, marker='*', color='y', s=35, alpha=0.5, edgecolor='k', label='Obs.')

    ax1.set_xlabel('sec.')
    ax1.set_ylabel('RMSE')
    ax1.set_xlim(0, 25)
    ax1.set_ylim(0, 10)
    ax1.set_title('Lorenz(1963) RMSE inflation +0.2', loc='left')

    plt.grid()
    plt.legend()
    plt.show()

class lorenz63_error_covariance_matrix:
  """
  一つのデータリストを誤差共分散行列の扱いを行うクラスメソッド
  """
  def __init__(self, path:str):
    self.err_data = self.read_error_csv(path)

  def read_error_csv(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")

  def error_heatmap(self, err_data:np.ndarray, timestep:int) -> None:
    fig, ax = plt.subplots()
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

if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  matrix_size  = 3
  obs_interval = 20
  mem          = 5000

  outdir    = './output/lorenz63/'
  da_method = 'KF'
  data_path = outdir + da_method + '.csv'
  err_path  = outdir + 'Error_matrix_' + da_method + '.csv'
  if da_method is 'EnKF':
    data_path = outdir +  da_method + '_' + str(mem) + 'm.csv'
    err_path  = outdir + 'Error_matrix_' + da_method + '_' + str(mem) + 'mem.csv'

  #---------------------------------------------------------- 
  # +++ class set
  # > lorenz63 cal. score
  score = lorenz63_score(data_path)
  # > prediction err covariance matrix
  e_matrix = lorenz63_error_covariance_matrix(err_path)
  
  #---------------------------------------------------------- 
  # +++ draw func.
  # >> 3D draw
  """
  time_before = time.time()
  laststep = 2500
  viewstep = int(laststep / 3) +1
  values = [ i_num for i_num in range(0, len(score.time_list[0:laststep]), 3)]
  with Pool(4) as pool:
    sleeps = list(tqdm(pool.imap(score.lorenz_3ddraw, values), total=viewstep))
  time_after = time.time()
  elapsed_time = time_after - time_before
  print(f'time: {elapsed_time} sec.')
  """

  # >> rmse draw
  num_sim_elem = 3
  rmse_sim = score.accuracy_rmse_func(score.true_list, score.sim_list, num_sim_elem)
  rmse_da  = score.accuracy_rmse_func(score.true_list, score.da_list, num_sim_elem)
  num_obs_elem = 2
  rmse_obs = score.accuracy_rmse_func(score.obs_true_list, score.obs_list, num_obs_elem)

  # Basic.
  score.lorenz_rmse_draw(rmse_sim, rmse_obs, rmse_da)

  # >> plus. dataset rmse
  """
  data_path_1 = outdir + 'EnKF_3m_0.2d0infla.csv'
  score_1 = lorenz63_score(data_path_1)
  data_path_2 = outdir + 'EnKF_10m_0.2d0infla.csv'
  score_2 = lorenz63_score(data_path_2)
  data_path_3 = outdir + 'EnKF_100m_0.2d0infla.csv'
  score_3 = lorenz63_score(data_path_3)
  data_path_4 = outdir + 'EnKF_1000m_0.2d0infla.csv'
  score_4 = lorenz63_score(data_path_4)
  data_path_5 = outdir + 'EnKF_5000m_0.2d0infla.csv'
  score_5 = lorenz63_score(data_path_5)

  rmse_da_1 = score.accuracy_rmse_func(score.true_list, score_1.da_list, num_sim_elem)
  rmse_da_2 = score.accuracy_rmse_func(score.true_list, score_2.da_list, num_sim_elem)
  rmse_da_3 = score.accuracy_rmse_func(score.true_list, score_3.da_list, num_sim_elem)
  rmse_da_4 = score.accuracy_rmse_func(score.true_list, score_4.da_list, num_sim_elem)
  rmse_da_5 = score.accuracy_rmse_func(score.true_list, score_5.da_list, num_sim_elem)

  score.lorenz_rmse_draw(rmse_sim, rmse_obs, rmse_da, rmse_da_1, rmse_da_2, rmse_da_3, rmse_da_4, rmse_da_5)
  """

  # >> error covariance matrix
  """
  for i_num in tqdm(range(0, int(len(score.time_list)/2), obs_interval)):
    matrix_data = e_matrix.err_data[i_num].reshape(matrix_size, matrix_size)
    e_matrix.error_heatmap(matrix_data, i_num)
  """