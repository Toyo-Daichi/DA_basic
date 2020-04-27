#-*- coding: utf-8 -*-
"""
Created on 2020.4.25
@author: Toyo_Daichi
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def df2_list(Dataframe: pd.core.frame.DataFrame, *index_wrd: tuple) -> list:
  index_list = []
  for i_index in index_wrd:
    index_list += Dataframe[i_index].tolist()
  return index_list

def read_Lorenz63_csv(path: str) -> list:
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

def read_error_csv(path: str) -> np.ndarray:
  return np.genfromtxt(path, delimiter=",")

def main_3ddraw(data: list):
  fig = plt.figure()
  ax = fig.gco(projection='3d')

def error_heatmap(err_data: np.ndarray):
  fig = plt.figure()
  cmap = sns.diverging_palette(220, 10, as_cmap=True)
  sns.heatmap(
    data=err_data, cmap=cmap, annot=True, fmt='g'
  )

  plt.show()

if __name__ == "__main__":
  #---------------------------------------------------------- 
  # +++ info. setting
  outdir    = './output/'
  da_method = 'EnKF'
  data_path = outdir + 'lorenz63_' + da_method + '.csv'
  err_path  = outdir + 'Error_matrix_lorenz63_' + da_method + '.csv'
  
  timescale   = 2500
  matrix_size = 3
  
  #---------------------------------------------------------- 
  # +++ lorenz63 cal. score
  #time_list, true_list, sim_list, da_list = read_Lorenz63_csv(data_path)
  
  # +++ prediction err covariance matrix
  err_list = read_error_csv(err_path)
  for i_num in len(err_list):
    print(i_num)
    err_data = err_list[i_num].reshape(matrix_size, matrix_size)
    error_heatmap(err_data)

