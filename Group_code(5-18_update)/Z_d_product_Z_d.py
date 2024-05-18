from multiprocessing import Value
from shared_code import abelian, get_sparse_d_tuple, compute_Sn, partition, make_all_products_from_multiset, print_vals
from group_operations import op_Z_d_Z_d
from math import comb, ceil, floor, factorial
from random import randint

# generates a set S \subset Z^d product Z_d 
# with |S| = k,
# where the number of elements with i in the 
# second component is distro_of_signs[i]
# thus, k = |S| = sum(distro_of_signs)
#
# ex: make_S(3, 6, (1, 3, 2)) will make a set 
# S = {(x1, 0), (x2, 1), (x3, 1), (x4, 1), (x5, 2), (x6, 2)}
# where the x_j are 3-tuples of sparse integers 
# notice there is 1 element with a 0 sign, 3 elements with a 1 sign
# and 2 elements with a 2 sign
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
def compute_s_n_with_sign_distribution(n, d, distribution):
  S = make_S(d = d, k = sum(distribution), distro_of_signs= distribution)
  s_n = compute_Sn(S, n, op_Z_d_Z_d)
  return len(s_n)

# formula for |S^n| when |S| = k, and there are t elements
# with a 1 in the 2nd component, k-t elements with a 0
# (i.e. no other signs)
# We can prove this formula and it matches with all simulations, 
# (see the test code)
def formula_when_all_one(n, k, t, d):
  # this is just the abelian case
  if t == 0:
    return comb(n + k -1, k - 1)
  
  if t == k:
    res = comb(ceil(n/d) + k - 1, k - 1) ** (n % d)
    res *= comb(floor(n/d) + k - 1, k - 1)**(d - (n % d))
    return res 

  res = comb(n + (k - t) - 1, k - t - 1)
  for i in range(1, n+1):
    cur_term = comb(ceil(i/d) + t - 1, t-1)**(i % d) 
    cur_term *= comb(floor(i/d) + t - 1, t-1)**(d - (i % d))
    num_forms = min(d, i+1)
    cur_term *= comb(n - i + num_forms * (k - t) - 1, num_forms * (k - t) - 1)
    res += cur_term
  return res

# We think the largest value of |S^n| is always achieved with just 0 and 1 
# in the second component (we DONT have a proof of this). 
# Thus we only need to check how many 1s to use, 
# i.e. try all possibly values of t from 0 to k in the above formula
def conjectured_gamma_k_n(k, n, d, noisy = False):
  best = -1
  achieved_at = -1
  for t in range(0, k+1):
    val = formula_when_all_one(n, k, t, d)
    if val > best:
      best = val
      achieved_at = t 
  
  if noisy:
    print("k = {}, n = {}, best at {} signs".format(k, n, achieved_at))
    
  return best


def test_t_vals(k, n, d, t_step):
  if k % t_step != 0:
    raise ValueError("t_step should evenly divde k")
  
  best = -1
  achieved_at = -1

  val_at_t = [-1] * (k//t_step + 1)

  # manually check up to 10 signs, we know the best is probably around here
  for t in range(0, 10):
    val = formula_when_all_one(n, k, t, d)
    if val > best:
      best = val
      achieved_at = t 

  # get percentages for higher values to get a sense
  for t in range(0, k+1, t_step):
    val = formula_when_all_one(n, k, t, d)
    val_at_t[t//t_step] = val
    if val > best:
      best = val
      achieved_at = t 
  
  for i in range(0, len(val_at_t)):
    print(f"with {i * t_step} negative signs, {val_at_t[i]/best * 100}% of the best")
  
  print_vals(n, k, d)
  print("k = {}, n = {}, best at {} signs".format(k, n, achieved_at))


test_t_vals(200, 200, 10, 1)

# takes about 10 seconds to run:
#conjectured_gamma_k_n(100, 5000, 5, noisy = True)


def get_approximation_formula(n, k, d):
  return factorial(n) / (factorial(n//d))**d * abelian(n, k)

def test_approx(n, k, d):
  if n % d != 0:
    raise ValueError("must have d | n")
  approx_num = get_approximation_formula(n, k, d)
  real_num = conjectured_gamma_k_n(k, n, d)
  print(real_num/approx_num)

def test_on_Z_d(n, k, d):
  def id_Z_d(d):
    return (tuple([0] * d), 0)
  # make a set with k-t elements have 0 in the second component, 
  # t elements having a 1 in the second component
  our_set = [(get_sparse_d_tuple(d), j % d) for j in range(0, k)]
  multi_sets = [list(p) for p in partition(n, k)]
  assert len(multi_sets) == abelian(n, k)
  total = 0
  
  print("This should be the sizeish of the multisets ")
  print(factorial(n) / (factorial(n//d))**d)

  biggest = -1
  for ms in multi_sets:
    all_prods = make_all_products_from_multiset(ms, our_set, id_Z_d(d), op_Z_d_Z_d)
    biggest = max(biggest, len(all_prods))

  print(f"biggest multiset was {biggest}")
  #assert total == formula_when_all_one(n, k, t, d)

#test_on_Z_d(10, 4, 2)

# Given an n, k, d try all possible distributions of signs, 
# return true if the best is achieved at some distribution that only 
# contains 0 or 1
# ex: only_zero_one_achives_max(5, 5, 3, True)
def only_zero_one_achives_max(n, k, d, noisy = False):

  def only_zero_one(distro):
    for i in range(2, d):
      if distro[i] != 0:
        return False
    return True 
  
  def contains_only_zero_one(list_of_distros):
    for distro in list_of_distros:
      if only_zero_one(distro):
        return True 
    return False

  biggest_s_n = -1
  achieved_at_distros = []
  all_sign_distros = partition(k, d)

  for distro in all_sign_distros:
    cur_val = compute_s_n_with_sign_distribution(n, d, distro)

    if noisy:
      print("With the distribution {} we get {}".format(distro, cur_val))

    if cur_val > biggest_s_n:
      biggest_s_n = cur_val
      achieved_at_distros = [distro]
    elif cur_val == biggest_s_n:
      achieved_at_distros.append(distro) 
  
  if noisy:
    print("We achieved the best values at these distros: {}".format(achieved_at_distros))

  return contains_only_zero_one(achieved_at_distros)

# run the above test on lots of values
# ex: show_we_only_need_zero_one(6, 4, 4, True)
def show_we_only_need_zero_one(n_upper, d_upper, k_upper, noisy = False):
  good = True
  for n in range (1, n_upper + 1):
    for k in range (1, k_upper + 1):
      for d in range(2, d_upper + 1):
        if noisy:
          print("n = {}, k = {}, d = {} we get".format(n, k, d), end = " ")
        if not only_zero_one_achives_max(n, k, d):
          if noisy:
            print("fail")
          good = False
        elif noisy:
          print("success")

  if not good:
    print("at least one time we needed a different sign distribution")
  else:
    print("success on every value")

  

  

# Conjecture: when k is fixed, n gets large we get an exponent of n^{d(k-1)}
# with a coefficient of 1/(d(k-1)!)
# I think we had a start of a proof with this. It matches with the simulation below,
# n has to be quite large for it to converge. 
# e.x. try exponent(1000, 5, 3) and exponent(10000, 5, 3)
def exponent(n, k, d):
  print("\'Actual\'")
  size_s_n = formula_when_all_one(n=n, k=k, d=d, t = 1)
  print(size_s_n/(n**(d * (k - 1))))
  
  print("Conjecture:")
  print(1/(factorial(d*(k - 1))))

exponent(10000, 15, 3)


def try_bound_lots():
  def try_some_bound(n, k, d):
    form_gamma = conjectured_gamma_k_n(k, n, d)
    bound = d**n * abelian(n, k)
    print(f"n = {n}, n = {k}, d = {d}")
    print(f"actual = {form_gamma}, bound = {bound}, actual/bound = {form_gamma/bound}")

  for n in range(90, 150, 10):
    for k in range(90, 150, 10):
      for d in range(2, 10, 2):
        try_some_bound(n, k, d)