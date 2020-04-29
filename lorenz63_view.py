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
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D

def df2_list(Dataframe:pd.core.frame.DataFrame, *index_wrd:tuple) -> list:
  index_list = []
  for i_index in index_wrd:
    index_list.append(Dataframe[i_index].tolist())
  return index_list

def read_Lorenz63_csv(path:str) -> list:
  """
  return list.shape
  list[0] = x score, list[1] = x score, list[2] = x score, 
  """
  df = pd.read_csv(path, sep=',')
  time_list = df2_list(df, ' timestep')
  true_list = df2_list(df, ' x_true', ' y_true', ' z_true')
  sim_list  = df2_list(df, ' x_sim', ' y_sim', ' z_sim')
  da_list   = df2_list(df, ' x_da', ' y_da', ' z_da')
  return  time_list, true_list, sim_list, da_list

def read_error_csv(path:str) -> np.ndarray:
  return np.genfromtxt(path, delimiter=",")

#def lorenz_3ddraw(true_data:list, sim_data:list, da_data:list, timestep:int):
def lorenz_3ddraw(timestep:int):
  print('start' + str(timestep))
  time_data, true_data, sim_data, da_data = read_Lorenz63_csv(data_path)

  fig = plt.figure()
  ax = fig.gca(projection='3d')

  ax.set_xlabel("X Axis"); ax.set_ylabel("Y Axis");  ax.set_zlabel("Z Axis")
  ax.plot(true_data[0], true_data[1], true_data[2], lw=0.5)

  if timestep != 'None':
    ax.scatter(true_data[0][timestep], true_data[1][timestep], true_data[2][timestep], marker='o', color='b', label='True')
    ax.scatter(sim_data[0][timestep], sim_data[1][timestep], sim_data[2][timestep], marker='o', color='r', label='Sim.')
    ax.scatter(da_data[0][timestep], da_data[1][timestep], da_data[2][timestep], marker='o', color='g', label='Data assim.')

    time = (timestep+1)*0.01
    ax.set_title('Lorenz(1963) - {:.2f} sec.'.format(time))
    plt.legend()
    plt.savefig('./figure/Lorenz_xyz_{:0>5}step.png'.format(timestep))
  plt.close('all')

  def lorenz_3ddraw_wrapper(args) -> None:
    return lorenz_3ddraw(*args)

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
  ax.set_title('Back ground Error Cov. - {:.1f} sec.'.format(time), loc='left')
  plt.savefig('./figure/error_heatmap_{:.1f}sec.png'.format(time))

  plt.close('all')

if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  matrix_size  = 3
  obs_interval = 20
  mem          = 10

  outdir    = './output/'
  da_method = 'EnKF'
  data_path = outdir + 'lorenz63_' + da_method + '.csv'
  err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '.csv'
  if (da_method == 'EnKF'):
    err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '_' + str(mem) + 'mem.csv'

  #---------------------------------------------------------- 
  # +++ basic info
  # > lorenz63 cal. score
  time_list, true_list, sim_list, da_list = read_Lorenz63_csv(data_path)
  # > prediction err covariance matrix
  err_list = read_error_csv(err_path)
  
  # draw func.
  values = [ i_num for i_num in range(0, len(time_list[0]))]
  with Pool(4) as pool:
    sleeps = [pool.map(lorenz_3ddraw, values)]
  
  """
  for i_num in tqdm(range(0, len(err_list), obs_interval)):
    err_data = err_list[i_num].reshape(matrix_size, matrix_size)
    error_heatmap(err_data, i_num)
  """
