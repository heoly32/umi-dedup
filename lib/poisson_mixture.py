import numpy as np
import math
import collections
from . import apportion_counts, umi_data

K_MAX = 10

# Store results into object
class BICResults:
  def __init__(self, bic, estimate, k):
    self.bic = bic
    self.estimate = estimate
    self.size = k

# Functions to compute likelihood
def dpois(x, mu):
  return x*np.log(mu) - mu - math.lgamma(x+1)

def mixture_dist(obs, param):
  n = obs.size
  k = param[0].size
  tmp = np.zeros( (n, k) )
  for j in range(n):
    tmp[j,...] = np.log(param[0])  + dpois(obs[j], param[1])
  if k == 1:
    return tmp
  else:
    return np.logaddexp.reduce(tmp, 1)

def likelihood(data, obs, param):
  return np.sum(data * mixture_dist(obs, param))

# BIC for selecting number of components
def BIC(data, obs, param):
    k = param[0].size
    return -2*likelihood(data, obs, param) + (2 * k - 1) * math.log(np.sum(data))

# EM-algorithm updates
def mixing_weights(obs, param):
    k = param[0].size
    n = obs.size
    output = np.zeros((n, k))
    for j in range(n):
      output[j,...] = (np.log(param[0]) + dpois(obs[j], param[1])) - np.logaddexp.reduce(np.log(param[0]) + dpois(obs[j], param[1]))
    return output

# def gradient(data, obs, param):
#   k = param[0].size
#   mixing_mat = mixing_weights(obs, param)
#   for j in range(k):
#     mixing_mat[...,j] += np.log(data)
#   gradient_prob = param[0] * np.logaddexp.reduce(mixing_mat, axis = 0)
#   gradient_theta = - np.logaddexp.reduce(mixing_mat, axis = 0)
#   mixing_mat = np.exp(mixing_mat)
#   for j in range(k):
#     mixing_mat[...,j] *= obs
#   gradient_theta += np.sum(mixing_mat, axis = 0)/param[1]
#   return np.concatenate((gradient_prob, gradient_theta))

def update_param(data, obs, param):
  k = param[0].size
  mixing_mat = np.exp(mixing_weights(obs, param))
  for j in range(k):
    mixing_mat[...,j] *= data
  next_prob =np.sum(mixing_mat, axis = 0)/np.sum(data)
  next_theta = 1/np.sum(mixing_mat, axis = 0)
  for j in range(k):
    mixing_mat[...,j] *= obs
  next_theta *= np.sum(mixing_mat, axis = 0)
  return (next_prob, next_theta)

# QN1 algorithm----
#We need to check if our proposed updates still lie in the parameter space
def in_param_space(current_param, param_step):
  k = current_param[0].size
  suggest_prob = current_param[0] + param_step[0:k]
  suggest_theta = current_param[1] + param_step[k:]
  return all(suggest_theta > 0) and all(suggest_prob > 0) and all(suggest_prob < 1)

#Define our updating step for the estimate of the inverse jacobian
def update_A(current_A, param_step, function_step):
  A_step = np.outer(param_step - np.dot(current_A, function_step), np.dot(param_step, current_A))
  A_step /=  np.dot(param_step, np.dot(current_A, function_step))
  return current_A + A_step

#We are conceptually looking for a zero of g_tilde
def g_tilde(data, obs, param):
  next_param = update_param(data, obs, param)
  next_step = (next_param[0] - param[0], next_param[1] - param[1])
  return np.concatenate(next_step)

# Fit algorithm to data----
def QN1_algorithm(data, obs, init_param):
    #parameter initialization
    next_param = init_param
    K = init_param[0].size
    next_A =  -np.identity(2*K)
    next_gtilde = g_tilde(data, obs, next_param)
    iter = 0
    if K == 1:
        next_prob = 1.0
        next_theta = float(np.sum(data * obs))/np.sum(data)
        next_param = (np.array(next_prob), np.array(next_theta))
    else:
        while True:
            iter += 1
            current_param = next_param
            current_gtilde = next_gtilde
            current_A = next_A
            #update parameter
            param_step = -np.dot(current_A, current_gtilde)
            #test if proposed parameter is in parameter space
            while not in_param_space(current_param, param_step):
                param_step /= 2
            #accept proposition
            next_param = (current_param[0] + param_step[0:K],
                                     current_param[1] + param_step[K:])
            #update the other parameters
            next_gtilde = g_tilde(data, obs, next_param)
            next_A = update_A(current_A, param_step,
                                             next_gtilde - current_gtilde)
            #testing if stopping rule is met
            current_lkhd = likelihood(data, obs, current_param)
            next_lkhd = likelihood(data, obs, next_param)
            if abs(current_lkhd - next_lkhd) < 10**-6:
              break
    # Compute BIC
    bic = BIC(data, obs, next_param)
    output = BICResults(bic, next_param, K)
    return output

def select_num_comp(data, obs):
  n = data.size
  bic_list = [QN1_algorithm(data, obs, (np.array(k*[1.0/k]), np.arange(1, k+1))) for k in range(1, min(K_MAX, n) + 1)]
  min_bic_result = min(bic_list, key = lambda p: p.bic)
  return min_bic_result

def dedup_cluster(umi_counts):
  naive_est = len(umi_counts.nonzero_keys())
  counter = collections.Counter(umi_counts.values())
  data = obs = []
  for count_item, data_item in counter.iteritems(): # PYTHON 2 ALERT
    obs.append(count_item)
    data.append(data_item)
  data = np.array(data)
  obs = np.array(obs)
  if data.size <= 2:
    est = naive_est
  else:
    min_bic_result = select_num_comp(data, obs)
    est = 0
    num_mol = np.arange(min_bic_result.size)
    mixing_mat = np.exp(mixing_weights(obs, min_bic_result.estimate))
    for i in range(data.size):
      est += np.dot(mixing_mat[i,...], num_mol) * obs[i]
    data_dedup = apportion_counts.apportion_umi_values(umi_counts, max(int(round(est)), naive_est))
    return umi_data.UmiValues(zip(umi_list, data_dedup))
