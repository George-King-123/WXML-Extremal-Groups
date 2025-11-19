from math import floor, ceil, factorial, sqrt
from shared_code import incl_range, format_large_num, choose, big_rand_num, compute_Sn
import group_operations

def compute_s_n_new_form(n, k, t): 
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
  

  lone_binomial = choose(n + k - t - 1, k - t - 1)
  first_sum = sum(N(k-t, ell) for ell in range(0, n) if n % 2 == ell % 2)
  
  def second_sum(): 
    def inner_sum(ell):
      return sum(R(t, q) * N(k-t, ell - q) for q in incl_range(1, ell))

    return sum(inner_sum(ell) for ell in incl_range(1, n) if n % 2 == ell % 2)

  return lone_binomial + first_sum + second_sum()

def compute_s_n_old_form(n, k, t):
  def I_0(t, k, n):
    return choose(n + (k - t) - 1, (k - t) - 1)

  def I_1(t, k, n):
    res = 0
    for q in range (1, min(k-t, n-1) + 1):
      res += pow(2, q)*choose(k-t, q)*choose(n-2, q-1)
    res *= t
    return res

  def I_2(t, k, n):
    res = 0
    for q in range (1, min(k-t, n-2) + 1):
      res += (pow(2, q) - 1) * choose(k-t, q) * choose(n-3, q-1) 
    # res *= t
    return res

  def R(t, k, n): 
    if t in {0, 1}: 
      return 0

    res = 0
    for i in range(2, n + 1):
      omega_value = 0
      sigma_value = 0
      for q in range (1, min(k - t, n - i) + 1):
        omega_value += (pow(2, q)) * choose(k-t, q) * choose(n-i -1, q-1)

      for m in range (1, min(t-1, ceil(i / 2)) + 1):
        sigma_value += choose(t, m) * choose(ceil(i / 2) - 1, m - 1) * choose((t-m) + floor(i / 2) - 1, (t - m) - 1)
      
      if i == n:
        omega_value = 1
      res += (omega_value * sigma_value)
    return res

  if t == 0:
    return choose(n + k - 1, k - 1)
  
  start = 1 if n % 2 == 0 else k
  cur = 2 if n % 2 == 0 else 3

  tot_I_0_contrib = 0 
  tot_I_1_contrib = 0
  tot_I_2_contrib = 0
  tot_R_contrib = 0

  while cur <= n: 
    tot_I_0_contrib += I_0(t, k, cur)
    tot_I_1_contrib += I_1(t, k, cur)
    tot_I_2_contrib += I_2(t, k, cur)
    tot_R_contrib += R(t, k, cur)
    cur += 2

  total = tot_I_0_contrib + tot_I_1_contrib + tot_I_2_contrib + tot_R_contrib + start
  
  return total

def compute_s_n_function_simulation(n, k, t): 
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

  S =  [-1] * t + [1] * (k-t)
  return len(run_simulation(S, n, k))

def compute_s_n_direct_simulation(n, k, t): 
  def get_S(t):
    first_components = [big_rand_num() for _ in range(k)]
    second_components = [-1] * t + [1] * (k-t)
    return zip(first_components, second_components)
  
  S = get_S(t)
  return len(compute_Sn(S=S, n=n, group_op=group_operations.op_d_inf))

def gamma_D_inf(n, k, noisy = False):
  biggest = -1
  achieved_at = 0
  for t in incl_range(0, k):
    val = compute_s_n_old_form(n, k, t)
    if val > biggest:
      biggest = val
      achieved_at = t

  if noisy:
    print(f"γ({n=}, {k=}) = {format_large_num(biggest)}, achieved with {achieved_at} signs")

  return biggest 


def growth_rate_one_sign_assumption(small_k, large_n):
  val = compute_s_n_old_form(n = large_n, k = small_k, t = 1)
  theta_bound = large_n**(small_k - 1)
  print(f"Actual ratio: {val/theta_bound}")

  coefficient = (2**(small_k-1))/factorial(small_k-1)
  print(f"Conjectured ratio: {coefficient}")

def gamma_n_n_exponent_two_sign_assumption(n):
  val = compute_s_n_old_form(n, n, t=2)
  experimental = val ** (1/n)
  print("experimental: {}".format(experimental))

  conjecture = 3 + 2 * sqrt(2)
  print("conjecture: {}".format(conjecture))

def check_new_form_matches_old(N = 30, k = 30, noisy = True): 
  for n in range(1, N + 1):
    for k in range(1, n + 1):
      for t in range(1, k+1): 
        new = compute_s_n_new_form(n, k, t)
        old = compute_s_n_old_form(n, k, t) 
        if noisy: 
          print(f"{new=} {old=}")
        assert new == old 

def check_limits_with_dif_t_vals(): 
  t = 14
  k = 14
  N = 1000
  s_n = compute_s_n_old_form(N, k, t)
  quotient = s_n / (N**(k-1))
  conj_limit = 2**(k - 2 * t + 1) * choose(2 * t - 2, t- 1) / factorial(k - 1)
  
  print(f"{conj_limit = }")
  print(f"{quotient = }")

def gamma_n_n_limit():
  N = 400
  alpha = .7
  t = round(N * alpha)
  s_n = compute_s_n_old_form(N, N, t)
  conj_limit = 3 + 2 * sqrt(2)  
  emprirical_limit = s_n**(1/N) 
  print(f"{conj_limit = }")
  print(f"{emprirical_limit = }")

def main():
  gamma_D_inf(n=20, k = 7, noisy=True)

if __name__ == "__main__":
  main()