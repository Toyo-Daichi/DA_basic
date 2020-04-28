#-*- coding: utf-8 -*-

import numpy as np
from sklearn.metrics import *

def mean_absolute_error(true:list, prediction:list) -> list:
  return mean_absolute_error(true, prediction)

def rmse(true:list, prediction:list) -> list:
  return np.sqrt(mean_squared_error(true, prediction))

def mean_squared_error(true:list, prediction:list) -> list:
  return mean_squared_error(true, prediction)

def correlation_coefficient(sig_xy:float, sig_x:float, sig_y:float) -> float:
  return sig_xy / (sig_x*sig_y)

