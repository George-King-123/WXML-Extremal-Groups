# run with python -m Z_d_product_Z_d

from shared_code import get_sparse_d_tuple, compute_Sn, incl_range, WC, format_large_num
from group_operations import op_Z_d_Z_d
from math import comb, ceil, floor, factorial
from tqdm import tqdm as loading_bar 
import time

# generates a set S \subset Z^d product Z_d with |S| = k,
# where the number of elements with i in the second component is distro_of_signs[i]
def make_S(k, d, distro_of_signs):
  if len(distro_of_signs) != d:
    raise ValueError("must say how many of each sign")
  
  if sum(distro_of_signs) != k:
    raise ValueError("total number of elements does not agree with the number of signs")
    
  S = set()
  for i in range(0, d):
    # add the right amount of elements of that sign
    for _ in range(0, distro_of_signs[i]):
      S.add((get_sparse_d_tuple(d), i))
  return S

# Returns |S^n| when S is constructed with the sign distrubtion as described
# in make_S
def compute_s_n_with_sign_distribution(n, distro_of_signs):
  S = make_S(d = len(distro_of_signs), k = sum(distro_of_signs), distro_of_signs= distro_of_signs)
  s_n = compute_Sn(S, n, op_Z_d_Z_d)
  return len(s_n)

def formula_when_all_zero_one_in_snd_comp(n, k, p, d):
  if p == 0:
    return comb(n + k - 1, k - 1)
  
  if p == k:
    res = comb(ceil(n/d) + k - 1, k - 1) ** (n % d)
    res *= comb(floor(n/d) + k - 1, k - 1)**(d - (n % d))
    return res 

  res = comb(n + (k - p) - 1, k - p - 1)
  for i in incl_range(1, n):
    cur_term = comb(ceil(i/d) + p - 1, p-1)**(i % d) 
    cur_term *= comb(floor(i/d) + p - 1, p-1)**(d - (i % d))
    num_forms = min(d, i+1)
    cur_term *= comb(n - i + num_forms * (k - p) - 1, num_forms * (k - p) - 1)
    res += cur_term
  return res


# test that the formula for just 0s and 1s in the second component
# works for all values of d,n,k up to the given upper bounds
# these have to be small because we explicitly compute |S^n| in order 
# to check the formula 
# ex: test_formula_all_one_many_vals(5, 5, 5)
def test_formula_all_one_many_vals(d_upper, n_upper, k_upper):
  print("testing that the Z_d formula works")
  def formula_works(n, k, p, d):
    S = make_S(d = d, k = k, distro_of_signs= [k - p, p] + ([0] * (d-2)))
    s_n = compute_Sn(S, n, op_Z_d_Z_d)
    return len(s_n) == formula_when_all_zero_one_in_snd_comp(n, k, p, d)
  

  for d in incl_range(2, d_upper):
    for n in incl_range (1, n_upper):
      for k in incl_range(1, k_upper):
        for p in incl_range(0, k):
          if not formula_works(n, k, p, d):
            print(f"FAIL {n=} {k=} {p=} {d=}")
            return
          
  print("success")

# test \lim_{n \to \infty} formula(t, d, k, n) / n^{d(k-1)}. 
# This didn't end up being important in the paper
# This call doesn't take too long: test_predicted_limit(big_number= 10 ** 6, k=10, p=4, d=7)
def test_predicted_limit(big_number, k, p, d):
    def analytic_limit(k, p, d):
      numerator = factorial(d * (p - 1))
      denominator = ((factorial(p - 1)) ** d) * (d ** (d * (p - 1))) * factorial(d * (k - 1))
      return numerator / denominator

    empirical_limit = formula_when_all_zero_one_in_snd_comp(n=big_number, k=k, p=p, d=d) / big_number ** (d * (k - 1))
    print(f"{empirical_limit=}")
    
    analytic_limit_value = analytic_limit(k, p, d)
    print(f"{analytic_limit_value=}")

def find_maximizing_t_val(n, k, d): 
    maximizers = []
    maximimum = -1 

    for p in incl_range(0, k): 
        cur_val = formula_when_all_zero_one_in_snd_comp(n=n, k=k, p=p, d=d)
        if cur_val > maximimum:
            maximizers = [p]
            maximimum = cur_val 
        elif cur_val == maximimum:
            maximizers.append(p) 
    
    return maximizers


def weak_composition_simulation(distro_of_signs, n):
  d = len(distro_of_signs)
  k = sum(distro_of_signs)

  S = [i for i in range(d) for _ in range(distro_of_signs[i])]

  id = (tuple([0] * (k * d)), 0)
  s_0 = { id }
  
  def get_s_n(n):
    if n == 0: 
      return s_0
    
    s_n_minus_one = get_s_n(n-1)
    result = set()

    for first_comp, shift in s_n_minus_one:
      for i in range(k): 
        first_comp_of_product = list(first_comp)
        first_comp_of_product[d * i + shift] += 1
        snd_comp_of_product = (shift + S[i]) % d

        result.add((tuple(first_comp_of_product), snd_comp_of_product))

    return result
  
  s_n = get_s_n(n) 

  return len(s_n)

def gamma_zd(n, k, d, noisy = False): 
  maximizing_distros = []
  best_value = -1

  all_sign_distros = WC(k, d) 
  for distro in all_sign_distros: 
    size_s_n = weak_composition_simulation(distro_of_signs=distro, n=n)

    if size_s_n > best_value: 
      best_value = size_s_n 
      maximizing_distros = [distro]
    elif size_s_n == best_value:
      maximizing_distros.append(distro)

  if noisy:
    print(f"{maximizing_distros=}")
    print(f"gamma_zd({n=}, {k=}) = {format_large_num(best_value)} for {d=}")

  return best_value, maximizing_distros

def test_conjecture(maxN, maxK, maxD): 
  def test_conjecture_one_val(n, k, d): 
    best_value, maximizing_distros = gamma_zd(n=n, k=k, d=d)
    return any(
      distro[0] + distro[1] == k for distro in maximizing_distros
    )

  for n in incl_range(1, maxN): 
    for k in incl_range(1, maxK): 
      for d in incl_range(2, maxD): 
        print(test_conjecture_one_val(n=n, k=k, d=d))

def test_conjecture2(maxN, maxK, maxD): 
  def test_conjecture_one_val(n, k, d): 
    all_distros = WC(k, d)
    for distro in all_distros: 
      if distro[0] + distro[1] != k: 
        val_of_distro = weak_composition_simulation(distro_of_signs=distro, n=n)

        modified_distro = tuple([distro[0], k-distro[0]] + [0] * (d-2))
        val_of_modified_distro = weak_composition_simulation(distro_of_signs=modified_distro, n=n)

        if val_of_modified_distro < val_of_distro and distro[0] != 0:
          print(f"conjecture FAILS on {distro}. {val_of_distro=} {val_of_modified_distro=}")

  for n in loading_bar(incl_range(1, maxN)): 
    for k in incl_range(1, maxK): 
      for d in incl_range(2, maxD): 
        test_conjecture_one_val(n=n, k=k, d=d)

def test_wc_sim(n, k, d): 
  all_sign_distros = WC(k, d) 
  for distro in loading_bar(all_sign_distros): 
    advanced_val = weak_composition_simulation(distro_of_signs=distro, n=n)
    simple_val = compute_s_n_with_sign_distribution(distro_of_signs=distro, n=n)

    if advanced_val != simple_val:
      print("FAIL")
      return
  
  print("all good")


def simulation_speed():
  N = 15
  k = 4
  d = 3

  def compute_all_vals(f):
    all_sign_distros = WC(k, d) 
    for distro in loading_bar(all_sign_distros): 
      f(distro_of_signs=distro, n=N)


  simulations = [weak_composition_simulation, compute_s_n_with_sign_distribution]
  for f in simulations: 
    print(f"Testing speed of {f.__name__}")

    start = time.perf_counter()
    compute_all_vals(f)
    end = time.perf_counter() 

    time_taken = end - start 
    print(f"{f.__name__} took {time_taken:.2f} seconds to compute value for {N=}")




def main():
  # test_formula_all_one_many_vals(6, 6, 5)
  # test_wc_sim(n=10, k=5, d=4)
  # simulation_speed()
  # weak_composition_simulation((5, 1, 0, 0), 7)
  # gamma_zd(n=5, k=15, d=4)
  test_conjecture2(maxN=6, maxK=6, maxD=5)
  pass


if __name__ == "__main__":
  main()