# run with python -m D_infinity.compute_s_n

from math import floor, ceil
from shared_code import incl_range, format_large_num, choose, big_rand_num, compute_Sn
import group_operations
from tqdm import tqdm as loading_bar
import time 

def new_formula(n, k, p): 
  if p == 0: 
    return choose(n + k - 1, k - 1)

  def R(x, y): 
    if y == 1: 
      return x 
    
    def term(r):
      return choose(x, r) * choose(y//2 - 1, r-1) * choose(ceil(y/2) + x - r - 1, x - r - 1)

    return sum(term(r) for r in incl_range(1, min(x-1, y//2)))
  
  def N(x, y):
    if y == 0: 
      return 1
    
    def term(s):
      return (2**s) * choose(x, s) * choose(y-1, s-1)
    
    return sum(term(s) for s in incl_range(1, min(x, y)))
  

  lone_binomial = choose(n + k - p - 1, k - p - 1)
  first_sum = sum(N(k-p, ell) for ell in range(0, n) if n % 2 == ell % 2)
  
  def second_sum(): 
    def inner_sum(ell):
      return sum(R(p, q) * N(k-p, ell - q) for q in incl_range(1, ell))

    return sum(inner_sum(ell) for ell in incl_range(1, n) if n % 2 == ell % 2)

  return lone_binomial + first_sum + second_sum()

def old_formula(n, k, p):
  def I_0(t, k, n):
    return choose(n + (k - p) - 1, (k - p) - 1)

  def I_1(p, k, n):
    res = 0
    for q in range (1, min(k-p, n-1) + 1):
      res += pow(2, q)*choose(k-p, q)*choose(n-2, q-1)
    res *= p
    return res

  def I_2(p, k, n):
    res = 0
    for q in range (1, min(k-p, n-2) + 1):
      res += (pow(2, q) - 1) * choose(k-p, q) * choose(n-3, q-1) 
    # res *= p
    return res

  def R(p, k, n): 
    if p in {0, 1}: 
      return 0

    res = 0
    for i in range(2, n + 1):
      omega_value = 0
      sigma_value = 0
      for q in range (1, min(k - p, n - i) + 1):
        omega_value += (pow(2, q)) * choose(k-p, q) * choose(n-i -1, q-1)

      for m in range (1, min(p-1, ceil(i / 2)) + 1):
        sigma_value += choose(p, m) * choose(ceil(i / 2) - 1, m - 1) * choose((p-m) + floor(i / 2) - 1, (p - m) - 1)
      
      if i == n:
        omega_value = 1
      res += (omega_value * sigma_value)
    return res

  if p == 0:
    return choose(n + k - 1, k - 1)
  
  start = 1 if n % 2 == 0 else k
  cur = 2 if n % 2 == 0 else 3

  total = start

  while cur <= n: 
    total += I_0(p, k, cur)
    total += I_1(p, k, cur)
    total += I_2(p, k, cur)
    total += R(p, k, cur)
    cur += 2
  
  return total

def function_simulation(n, k, p): 
  def run_simulation(S, n, k):

    def compute_one_more(cur_tuple, last_fcn_value, S):
      # either +1 or -1
      cur_sign = cur_tuple[1]

      # cur_sum is a list of length k indicating, where 
      # cur_sum[i] is the (signed) number of times the ith 
      # elt of S appears in in the sum
      cur_sum = list(cur_tuple[0])
      cur_sum[last_fcn_value] += cur_sign
      cur_sign *= S[last_fcn_value]

      return (tuple(cur_sum), cur_sign)

    if n == 0:
      return {(tuple([0] * k), 1)}
    len_n_minus_one_prods = run_simulation(S, n-1, k)
    len_n_prods = set()
    for prod in len_n_minus_one_prods:
      for i in range(0, k):
        len_n_prods.add(compute_one_more(prod, i, S))
    # last fcn_value is f(n) for our current value of n

    return len_n_prods

  S =  [-1] * p + [1] * (k-p)
  return len(run_simulation(S, n, k))

def direct_simulation(n, k, p): 
  def get_S(p):
    first_components = [big_rand_num() for _ in range(k)]
    second_components = [-1] * p + [1] * (k-p)
    return set(zip(first_components, second_components))
  
  S = get_S(p)
  return len(compute_Sn(S=S, n=n, group_op=group_operations.op_d_inf))

# check all p <= k <= n, up to n = k = p = N
def match_on_all_vals(f1_nkp, f2_nkp, N, print_all = False): 
  all_match = True
  for n in loading_bar(incl_range(1, N)):
    for k in incl_range(1, n): 
      for p in incl_range(0, k): 
        f1_val = f1_nkp(n=n, k=k, p=p)
        f2_val = f2_nkp(n=n, k=k, p=p)
        if f1_val != f2_val:
          all_match = False
          print(f"FAIL on {n=}, {k=}, {p=}: f1 = {format_large_num(f1_val)}, f2 = {format_large_num(f2_val)}")
        
        if print_all:
          print(f"{n=} {k=} {p=}, f1 = {format_large_num(f1_val)} f2 = {format_large_num(f2_val)}")
  
  return all_match

all_ways_of_computing_s_n = [new_formula, old_formula, function_simulation, direct_simulation]

def check_all_ways_of_computing_s_n_same():
  # simulations are quite slow, don't set N much higher than this
  N = 8

  print(f"Checking that all {len(all_ways_of_computing_s_n)} ways of computing S^n match,"
      + f"for all p <= k <= n, up to n = {N}")
  
  all_consecutive_pairs = zip(all_ways_of_computing_s_n, all_ways_of_computing_s_n[1:])

  for f1, f2 in all_consecutive_pairs:
    print(f"Checking {f1.__name__} against {f2.__name__}")

    if match_on_all_vals(f1, f2, N=N): 
      print(f"{f1.__name__} matches {f2.__name__}")
    else: 
      print(f"FAIL: {f1.__name__} DOES NOT MATCH {f2.__name__}")

    separator_for_readability = ("=" * 20)
    print(separator_for_readability) 

def check_formulas_same(): 
  N = 40
  print(f"Checking that formulas match on all values, for all t <= k <= n, up to n = {N}")
  match_on_all_vals(f1_nkp=new_formula, f2_nkp=old_formula, N=N)

# used to figure out which formula was faster, and thus should be the canonical formula
# For N = 50, old takes 21 seconds and new takes 27 seconds 
def test_formula_speed():
  N = 50

  def compute_all_vals(f):
    for n in loading_bar(incl_range(1, N)):
      for k in incl_range(1, n):
        for p in incl_range(0, k):
          f(n=n, k=k, p=p)

  formulas = [old_formula, new_formula]
  for f in formulas: 
    print(f"Testing speed of {f.__name__}")

    start = time.perf_counter()
    compute_all_vals(f)
    end = time.perf_counter() 

    time_taken = end - start 
    print(f"{f.__name__} took {time_taken:.2f} seconds to compute value for {N=}")

FAST_FORMULA = old_formula

def gamma_D_inf(n, k):
  return max(FAST_FORMULA(n, k, p) for p in incl_range(0, k))

def find_maximizing_p_values(n, k):
    maximizers = []
    maximimum = -1 

    for p in incl_range(0, k): 
        cur_val = FAST_FORMULA(n=n, k=k, p=p)
        if cur_val > maximimum:
            maximizers = [p]
            maximimum = cur_val 
        elif cur_val == maximimum:
            maximizers.append(p) 
    
    return maximizers
  
def main():
  # check_all_ways_of_computing_s_n_same()
  check_formulas_same()
  # test_formula_speed()


if __name__ == "__main__":
  main()




