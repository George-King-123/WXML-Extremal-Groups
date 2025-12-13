# run with python -m Z_d_wreath_product_S_d

from math import factorial
import random 
from group_operations import op_wreath_prod
from shared_code import compute_Sn, get_sparse_d_tuple, WC, format_large_num
import Z_d_product_Z_d
from tqdm import tqdm as loading_bar

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


  iter_obj = loading_bar(range(num_iters)) if noisy else range(num_iters)

  for _ in iter_obj:
    S = get_rand_sparse_set(k = k, d = d)
    val = len(compute_Sn(S = S, n = n, group_op=op_wreath_prod))
    if val > best:
      best = val
      best_set = S 
  
  if noisy:
    print(f"Best value after trying {num_iters} sets: {format_large_num(best)}")
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
    # formula is relatively well known
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

# Since Z_d is a special case of S_d, we can test the S_d code by making sure a set 
# that only has "shift by one" permutations agrees with the formula we came up 
# with for Z_d
# using the given n, k, d, and t many components with 1 in the 2nd component
# with k - t elements having 0 in the 2nd component
def restrict_to_Z_d(n, k, d, p):
  def get_s_d_set():
    # shift_by_one \in Z_d
    shift_by_one = tuple(list(range(1, d)) + [0])
    identity = tuple(range(0, d))
    
    return set(
     [(get_sparse_d_tuple(d), shift_by_one) for _ in range(p)] + 
     [(get_sparse_d_tuple(d), identity) for _ in range(k-p)]  
    )
    
  val_from_s_d_sim = len(compute_Sn(S = get_s_d_set(), n = n, group_op = op_wreath_prod))
  val_from_z_d_formula = Z_d_product_Z_d.formula_when_all_zero_one_in_snd_comp(n=n, k=k, p=p, d=d)

  return val_from_s_d_sim == val_from_z_d_formula

# ex: test_a_bunch_of_restrictions(6, 6, 3)
def test_a_bunch_of_restrictions(n_upper, k_upper, d_upper):
  print("Testing S_d code by seeing that it lines up with Z_d")
  for n in range(1, n_upper + 1):
    for k in range(1, k_upper + 1):
      for d in range(1, d_upper + 1):
        for p in range(0, k + 1):
          if not restrict_to_Z_d(n=n, k = k, d= d, p = p):
            print("fail")

  print("success")


def main():
  #  experiment_with_d_3_bounds()
  # print_s_n_paper_set()
  see_how_tight_bounds_are()
  # test_a_bunch_of_restrictions(6, 6, 3)

if __name__ == "__main__": 
    main()
