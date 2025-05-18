from math import factorial, comb
import random 
from group_operations import op_wreath_prod
from shared_code import compute_Sn, get_sparse_d_tuple
import itertools

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
  assert len(res) == factorial(d)
  return res 


# returns S \subset Z^d wreathproduct S_d where |S| = k
# the numbers in the first component are sparse and the 
# permutations in the second component are random
def get_rand_sparse_set(k, d):
  s_d = list(make_s_d(d))  
  S = set()
  
  # add something to the set k times 
  for _ in range(0, k):
    rand_s_d_elt = s_d[random.randrange(factorial(d))]
    S.add((get_sparse_d_tuple(d), rand_s_d_elt))
  
  return S

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
    print(f"Best value after trying {num_iters} sets: {best:_}")
    print("The best set had these 2nd components:")
    print_second_components(best_set)

  return best

# just print the second component, we know the first component is d sparse integers
def print_second_components(S):
  print([elt[1] for elt in S])

# return the set mentioned in the paper, with identity in the second component except for (12) and (123...d)
def get_paper_set(k, d):
    perm_id = tuple([i for i in range(0, d)])
    perm_12 = tuple([1, 0] + list(range(2, d))) 
    perm_cyc = tuple([i + 1 for i in range(0, d - 1)] + [0])
    
    S = set()
    for i in range(0, k - 2):
        S.add((get_sparse_d_tuple(d), perm_id))
    S.add((get_sparse_d_tuple(d), perm_12))
    S.add((get_sparse_d_tuple(d), perm_cyc))
    return S

def compute_Sn_paper_set(n, k, d):
    S = get_paper_set(k, d)
    s_n = compute_Sn(S, n, op_wreath_prod)
    return len(s_n)

# h_n(r) is the number of n \times n magic squares with sum r, as defined in the overleaf
def h_3(r):
    return (r+1) * (r+2) * (r**2 + 3*r + 4) // 8

# f_d(k, n) = \sum_{x \in \WC(n, k)} \prod_{i=1}^{k} h_d(x_i) 
def f_3(k, n):
    WC_n_k = WC(n, k)
    tot = 0
    for x in WC_n_k:
        prod = 1 
        for i in range(0, k):
            prod *= h_3(x[i])
        tot += prod
    
    return tot

def WC_generator(n, k):
  if k <= 0 or n < 0:
     raise ValueError()
  
  for bars in itertools.combinations(range(n + k - 1), k - 1):
      result = []
      prev = -1
      for b in bars:
          result.append(b - prev - 1)
          prev = b
      result.append(n + k - 1 - prev - 1)
      yield result

def WC(n, k):    
  wc = list(WC_generator(n, k)) 
  assert len(wc) == comb(n + k - 1, k - 1) 
  assert all(sum(x) == n for x in wc)
  assert all(len(x) == k for x in wc)
  return wc

def lower_bound_3(k, n):
  tot = 0
  C_3 = 6
  for i in range(C_3, n):
     tot += (f_3(k-2, n-i) * (i - C_3 + 1))
  return tot

def see_how_tight_bounds_are():
  k=3
  n=15
  d=3
  print(f"Upper bound: f_3(k, n) = {f_3(k, n):_}")
  print(f"Lower bound: {lower_bound_3(k,n):_}")
  print(f"S^n paper set: {compute_Sn_paper_set(n=n, k=k, d=d):_}")
  compute_approx_gamma_k_n(k = k, n = n, d = d, num_iters = 100, noisy = True)

def experiment_with_d_3_bounds(): 
  
  def upper_bound(k, n):
    return f_3(k, n)
  
  K = 6
  N = 70

  print(f"k = {K}, n = {N}")

  upper = upper_bound(K, N)
  lower = lower_bound_3(K, N)

  print(f"Quotient of bounds: {(upper/lower):_}")

def print_s_n_paper_set():
  k = 4
  n = 15
  d = 3

  print(f"Upper bound: f_3(k, n) = {f_3(k, n):_}")
  print(f"Lower bound = {lower_bound_3(k,n):_}")
  print(f"S^n paper set: {compute_Sn_paper_set(n=n, k=k, d=d):_}")

def main():
  #  experiment_with_d_3_bounds()
  # print_s_n_paper_set()
  # see_how_tight_bounds_are()
  
  k = 4
  for n in range(10, 2000, 50):
    print(lower_bound_3(k, n) / (n ** (5 * k - 9)))

if __name__ == "__main__": 
    main()
