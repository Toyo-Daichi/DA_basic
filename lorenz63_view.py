#-*- coding: utf-8 -*-
"""
Created on 2020.4.25
@author: Toyo_Daichi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def df2_list(Dataframe: pd.core.frame.DataFrame, *index_wrd: tuple) -> list:
  index_list = []
  for i_index in index_wrd:
    index_list += Dataframe[i_index].tolist()
  return index_list

def read_csv(path: str) -> list:
  """
  return list.shape
  list[0] = x score, list[1] = x score, list[2] = x score, 
  """
  df = pd.read_csv(path)
  time_list = df2_list(df, ' timestep')
  true_list = df2_list(df, ' x_true', ' y_true', ' z_true')
  sim_list  = df2_list(df, ' x_sim', ' y_sim', ' z_sim')
  da_list   = df2_list(df, ' x_da', ' y_da', ' z_da')
  return  time_list, true_list, sim_list, da_list

def main_draw(data):
  fig = plt.figure()
  ax = fig.gco(projection='3d')

if __name__ == "__main__":
  # info. setting
  outdir    = './output/'
  da_method = 'EnKF'
  data_path = outdir + 'lorenz63_' + da_method + '.csv'

  time_list, true_list, sim_list, da_list = read_csv(data_path)