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
    self.df = pd.read_csv(path, sep=',')
    # sim. list
    self.time_list = self.df2_list(self.df, ' timestep')[0]
    self.true_list = self.df2_list(self.df, ' x_true', ' y_true', ' z_true')
    self.sim_list  = self.df2_list(self.df, ' x_sim', ' y_sim', ' z_sim')
    self.da_list   = self.df2_list(self.df, ' x_da', ' y_da', ' z_da')
    # obs. list
    undef = -999.e0
    obs = self.df.loc[ :, [' timestep', ' x_true', ' y_true', ' x_obs', ' y_obs']]
    obs = obs[obs[' x_obs'] != undef]
    self.obs_time_list  = self.df2_list(obs, ' timestep')[0]
    self.obs_true_list = self.df2_list(obs, ' x_true', ' y_true')
    self.obs_list       = self.df2_list(obs, ' x_obs', ' y_obs')
    """
    list.shape
    list[0] = x score, list[1] = y score, list[2] = z score, 
    """
    print('Call ..... Constract of lorenz63_score')
  
  def df2_list(self, Dataframe:pd.core.frame.DataFrame, *index_wrd:tuple) -> list:
    index_list = []
    for i_index in index_wrd:
      index_list.append(Dataframe[i_index].tolist())
    return index_list

  def accuracy_rmse_func(self, true_list:list, asses_list:list, num_elem:int) -> list:
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

  def lorenz_rmse_draw(self, rmse_sim:list, rmse_da:list, rmse_obs:list, da_method) -> None:
    fig = plt.figure()
    ax1 = fig.subplots()
    
    sns.set_style('whitegrid')
    ax1.plot(self.time_list, rmse_da, ls="--", color='r', label='Data assim.')
    ax1.scatter(self.obs_time_list, rmse_obs, marker='o', color='g', s=20, alpha=0.5, edgecolor='k', label='Obs.')

    ax1.set_xlabel('sec.')
    ax1.set_ylabel('RMSE')
    ax1.set_ylim(0, 1.6)
    ax1.set_xlim(0, 25)
    ax1.set_title('Lorenz(1963) RMSE, Data assim. method: ' + da_method, loc='left')

    plt.grid()
    plt.legend()
    plt.show()

class lorenz63_error_covariance_matrix:
  def __init__(self, path:str):
    self.err_data = self.read_error_csv(path)

  def read_error_csv(self, path:str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",")

  def error_heatmap(self, err_data:np.ndarray, timestep:int):
    fig, ax = plt.subplots()
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.set(style='white')
    sns.heatmap(
      data=err_data, cmap=cmap, annot=True, fmt='.6f',
      vmax=0.1, vmin=-0.1, center=0,
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

  outdir    = './output/'
  da_method = 'EnKF'
  data_path = outdir + 'lorenz63_' + da_method + '_little_diff.csv'
  err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '.csv'
  if da_method is 'EnKF':
    data_path = outdir + 'lorenz63_' + da_method + '_' + str(mem) + 'mem_little_diff.csv'
    err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '_' + str(mem) + 'mem.csv'

  #---------------------------------------------------------- 
  # +++ class set
  # > lorenz63 cal. score
  score = lorenz63_score(data_path)
  # > prediction err covariance matrix
  e_matrix = lorenz63_error_covariance_matrix(err_path)
  
  #---------------------------------------------------------- 
  # +++ draw func.
  # > 3D draw
  """
  time_before = time.time()
  laststep = 2500
  viewstep = int(laststep / 3) +1
  values = [ i_num for i_num in range(0, len(score.time_list[0][0:laststep]), 3)]
  with Pool(4) as pool:
    sleeps = list(tqdm(pool.imap(score.lorenz_3ddraw, values), total=viewstep))
  time_after = time.time()
  elapsed_time = time_after - time_before
  print(f'time: {elapsed_time} sec.')
  
  # > error covariance matrix
  for i_num in tqdm(range(0, len(err_list), obs_interval)):
    matrix_data = e_matrix.err_data[i_num].reshape(matrix_size, matrix_size)
    e_matrix.error_heatmap(matrix_data, i_num)
  """
  
  # > rmse draw
  num_sim_elem = 3
  rmse_sim = score.accuracy_rmse_func(score.true_list, score.sim_list, num_sim_elem)
  rmse_da  = score.accuracy_rmse_func(score.true_list, score.da_list, num_sim_elem)
  num_obs_elem = 2
  rmse_obs = score.accuracy_rmse_func(score.obs_true_list, score.obs_list, num_obs_elem)

  score.lorenz_rmse_draw(rmse_sim, rmse_da, rmse_obs, da_method)