#-*- coding: utf-8 -*-
"""
Created on 2020.4.25
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './lib'))

import cal_statics
import numpy as np
import pandas as pd
import seaborn as sns
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D

class lorenz63_score:
  def __init__(self, path):
    self.df = pd.read_csv(path, sep=',')
    self.time_list = self.df2_list(self.df, ' timestep')
    self.true_list = self.df2_list(self.df, ' x_true', ' y_true', ' z_true')
    self.sim_list  = self.df2_list(self.df, ' x_sim', ' y_sim', ' z_sim')
    self.da_list   = self.df2_list(self.df, ' x_da', ' y_da', ' z_da')
    """
    list.shape
    list[0] = x score, list[1] = x score, list[2] = x score, 
    """
    print('Call ..... Constract of lorenz63_score')
  
  def df2_list(self, Dataframe:pd.core.frame.DataFrame, *index_wrd:tuple) -> list:
    index_list = []
    for i_index in index_wrd:
      index_list.append(Dataframe[i_index].tolist())
    return index_list

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


def read_error_csv(path:str) -> np.ndarray:
  return np.genfromtxt(path, delimiter=",")

def error_heatmap(err_data:np.ndarray, timestep:int):
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
  mem          = 10

  outdir    = './output/'
  da_method = 'KF'
  data_path = outdir + 'lorenz63_' + da_method + '_little_diff.csv'
  err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '.csv'
  if (da_method == 'EnKF'):
    err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '_' + str(mem) + 'mem.csv'

  #---------------------------------------------------------- 
  # +++ basic info
  # > lorenz63 cal. score
  score = lorenz63_score(data_path)
  # > prediction err covariance matrix
  err_list = read_error_csv(err_path)
  
  # draw func.
  time_before = time.time()
  laststep = 2500
  viewstep = int(2500 / 3) +1
  values = [ i_num for i_num in range(0, len(score.time_list[0][0:laststep]), 3)]
  with Pool(4) as pool:
    sleeps = list(tqdm(pool.imap(score.lorenz_3ddraw, values), total=viewstep))
  time_after = time.time()
  elapsed_time = time_after - time_before
  print(f'time:{elapsed_time}s')
  
  """
  for i_num in tqdm(range(0, len(err_list), obs_interval)):
    err_data = err_list[i_num].reshape(matrix_size, matrix_size)
    error_heatmap(err_data, i_num)
  """
