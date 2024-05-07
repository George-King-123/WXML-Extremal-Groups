from math import factorial
from random import randint
from group_operations import op_wreath_prod
from shared_code import compute_Sn, get_sparse_d_tuple
from Z_d_product_Z_d import conjectured_gamma_k_n

# This file is less well documented than the others, be advised


# returns a set of all tuples of length d having elements in
# {0, 1, ..., d-1} representing bijections
def make_s_d(d:int):
  def select(cur, all_elts:set, available_nums:set):
    if len(available_nums) == 0:
      all_elts.add(tuple(cur))
      return
    
    for i in available_nums:
      cur.append(i)
      select(cur, all_elts, available_nums.difference({i}))
      cur.remove(i)
  
  res = set()
  select([], res, {i for i in range (0, d)})
  return res 


# returns S \subset Z^d wreathproduct S_d where |S| = k
# the numbers in the first component are sparse and the 
# permutations in the second component are random
def get_rand_sparse_set(k, d):
  s_d = list(make_s_d(d))

  if len(s_d) != factorial(d):
    raise ValueError("Something went horribly wrong")
  
  S = set()
  
  # add something to the set k times 
  for _ in range(0, k):
    rand_s_d_elt = s_d[randint(0, factorial(d) - 1)]
    S.add((get_sparse_d_tuple(d), rand_s_d_elt))
  
  return S


def test(k, n, d):
  S = get_rand_sparse_set(k=k, d=d)
  print(S)

  s_n = compute_Sn(S=S, n=n, group_op=op_wreath_prod)
  print(len(s_n))


# ex: compare_many_vals(d = 3, fixed_k = 4, n_lower = 10, n_upper = 25)
def compare_many_vals(d, n_lower, n_upper, fixed_k):
  S = get_rand_sparse_set(k=fixed_k, d=d)

  print_second_components(S)

  print("k = {}".format(fixed_k))
  for n in range (n_lower, n_upper + 1):
    z_d_val = conjectured_gamma_k_n(k=fixed_k, n=n, d=d)
    s_d_val = len(compute_Sn(S=S, n=n, group_op=op_wreath_prod))
    print("for n = {}, Z_d was {} and S_d was {}, ratio S_d/Z_d is {}".format(n, z_d_val, s_d_val, s_d_val/z_d_val))



# try num_iters different sets S with cardinality k, 
# compute |S^n| for each and take the max
# ex: compute_approx_gamma_k_n(k = 3, n = 8, d = 3, num_iters = 20, noisy = True)
def compute_approx_gamma_k_n(k, n, d, num_iters, noisy = False):
  best = -1
  best_set = set()

  for _ in range(0, num_iters):
    S = get_rand_sparse_set(k = k, d = d)
    val = len(compute_Sn(S = S, n = n, group_op=op_wreath_prod))
    if val > best:
      best = val
      best_set = S 
  
  if noisy:
    print("for n = {}, k = {}, d = {}".format(n, k, d) + 
          " trying {} different sets we get a value ".format(num_iters) +
           "of {}".format(best))
    print("The best set had these 2nd components:")
    print_second_components(best_set)

  return best


# just print the second component, we know the first component is d sparse integers
def print_second_components(S):
  print([elt[1] for elt in S])

def find_a_good_set(d, k, n, num_iters):
  best_set = set()
  best_val = 0
  prs = []

  value_set = set()

  for _ in range(0, num_iters):
    S = get_rand_sparse_set(k=k, d=d)

    num_distinct_second_components = len({elt[1] for elt in S})
    cur_val = len(compute_Sn(S = S, n = n, group_op= op_wreath_prod)) 
    
    value_set.add(cur_val)
    prs.append((cur_val, num_distinct_second_components))

    #print("With {} 2nd components we got {}".format(num_distinct_second_components, cur_val))
    if cur_val > best_val:
      best_set = S 
      best_val = cur_val 

  print("we got {} different values".format(len(value_set)))
  print(best_val)
  #print(sorted(prs))
  print_second_components(best_set)



#find_a_good_set(d = 3, k = 3, n = 15, num_iters=10)


compute_approx_gamma_k_n(k = 3, n = 8, d = 3, num_iters = 20, noisy = True)

def something(n_l, n_h, k_l, k_h, d, num_iters):
  for k in range (k_l, k_h + 1):
    for n in range(n_l, n_h + 1):
      s_d = compute_approx_gamma_k_n(k=k, n=n, d=d, num_iters=num_iters)
      z_d = conjectured_gamma_k_n(k=k, n=n, d=d)

      print("n = {}, k = {}, s_d/n^(d(k-1)) = {}".format(n, k, z_d/(n**(d*(k-1)))))

#something(5, 10, 2, 6, 3, 15)
