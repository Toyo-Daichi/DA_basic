# -*- coding:utf-8 -*-
# Reference: 
# https://satomacoto.blogspot.com/2012/11/python.html

import numpy as np

def model(x:np.array,u:np.array):
  """Transition model
  Args:
    x (np.array): var
    u (np.array): parm
  Return:
    x_transition(np.array): time integral
  Note:
    x(t=i+1) = f(x(t=i),u)*dt
             = (Ax(t=1) + Bu + noise)*dt, noise ~ N(0,2I)
    ** 
    Matrix size
    [2,1]    = [2,2]*[2,1] + [2,2]*[2,1] + [2,1] 
    **
    For simplicity, assume that dt = 1.0s.
  """
  # parm matrix
  var, err = 2, 2
  A = np.eye(2)
  B = np.eye(2)
  
  # make normalized random number
  noise_mean = np.zeros(var)
  noise_cov  = np.eye(var)*err
  noise = np.random.multivariate_normal(noise_mean,noise_cov)
  
  # calculate transition
  x_transition = A@x + B@u + noise
  
  return x_transition

def observation(x:np.array, err):
  """Observation model
  Args:
    x(np.array): var 
  Return:
    x_obs(np.array): observation
  """
  # observation point
  var = 4
  point_obs = (
    [ 0.,  0.], [10.,  0.],
    [ 0., 10.], [10., 10.]
  )

  # make normalized random number
  noise_mean = np.zeros(var)
  noise_cov  = np.eye(var)*err
  noise = np.random.multivariate_normal(noise_mean,noise_cov)

  # observation vactor
  x_obs = np.array([
    np.linalg.norm(x-point_obs[0]), np.linalg.norm(x-point_obs[1]),
    np.linalg.norm(x-point_obs[2]), np.linalg.norm(x-point_obs[3])
  ]) + noise

  return x_obs

def likelihood(x,y):
  """Likelihood function """
  likelihood = np.exp(-((y-observation(x,0)))@(y-observation(x,0))/4)
  return likelihood

def importance_sampling(weight,x_fcst,x_obs):
  weight_update = weight*likelihood(x_fcst,x_obs)
  return weight_update

def re_sampling(x,weight):
  """Re-sampling """
  x_resampled = np.zeros_like(x)
  w_resampled = np.zeros_like(weight)
  sample = np.random.multinomial(2,weight)
  for i, n in enumerate(sample):
    for j in range(n):
      x_resampled[i] = x[i]
  w_resampled[:] = 1.0
  return x_resampled, w_resampled

def particle_filter(x_step,x_obs,weight,particle,parm):
  """ DA(particle filter) """
  for i_particle in range(particle):
    x_step[i_particle,:] = model(x_step[i_particle], parm)
    weight[i_particle] = importance_sampling(weight[i_particle],x_step,x_obs)
  
  # normalization
  weight_sum = np.sum(weight)
  for i_particle in range(particle):
    weight[i_particle] /= weight_sum
  
  x_resampled, weight_resampled = re_sampling(x_step,weight)

  return x_resampled,weight_resampled
   
if __name__ == '__main__':
  # set initial parm.
  step     = 10
  particle = 10
  err      = 4.0
  #
  x_init   = np.array([0.,0.]) # <= two variables 
  parm     = np.array([4.,4.])

  """ make true & obs """
  x_true = np.zeros((step,2))
  x_obs  = np.zeros((step,4))
  x = x_init
  for it in range(step):
    x = model(x,parm)
    y = observation(x,err)
    x_true[it,:], x_obs[it,:] = x, y

  """ forecast """
  x_particle = np.zeros((step,particle,2))
  x_anl = np.zeros((step,2))
  weight = np.zeros((step,particle))

  # make initial score & weight
  for i_particle in range(particle):
    x_particle[0,i_particle,:] = np.random.randn(2)*20
    weight[0,i_particle] = 1.0

  """ DA(particle filter) """
  for it in range(step-1):
    x_anl[it,0], x_anl[it,1] = np.mean(x_particle[it,:,0]), np.mean(x_particle[it,:,1])
    print('TRUE:{:.3f}, FCST:{:.3f}, RMSE:{:.3f}'.format(x_true[it,0],x_anl[it,0],((x_true[it,0]-x_anl[it,0])**2)**(0.5)))

    # main
    x_particle[it+1,:,:], weight[it+1] = particle_filter(x_particle[it],x_obs[it],weight[it],particle,parm)
