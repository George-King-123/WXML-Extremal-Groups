from math import factorial, comb, log
from random import randint
from re import I
from typing import OrderedDict
from group_operations import op_s_d, op_wreath_prod
from shared_code import abelian, compute_Sn, get_sparse_d_tuple, partition, make_all_products_from_multiset, print_vals
from Z_d_product_Z_d import conjectured_gamma_k_n


# This file is less well documented than the others, be advised

def get_id_wreath(d):
  return (tuple([0] * d), tuple([i for i in range(0, d)]))

# returns a set of all tuples of length d having elements in
# {0, 1, ..., d-1} representing bijections
def make_s_d(d):
  def select(cur, all_elts, available_nums):
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

def generate_all_transpositions_in_s_d(d):
  res = set()

  # a transposition switches i and j 
  for i in range(0, d):
    for j in range(i + 1, d):
      sigma = [i for i in range(0, d)]
      sigma[i] = j
      sigma[j] = i 
      res.add(tuple(sigma))
  
  assert len(res) == comb(d, 2)
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

# just print the second component, we know the first component is d sparse integers
def print_second_components(S):
  print([elt[1] for elt in S])

# just picks one random set, not that good 
# ex: compare_many_vals(d = 3, fixed_k = 4, n_lower = 10, n_upper = 25)
def compare_many_vals(d, n_lower, n_upper, fixed_k):
  S = get_rand_sparse_set(k=fixed_k, d=d)

  print_second_components(S)

  print("k = {}".format(fixed_k))
  for n in range (n_lower, n_upper + 1):
    z_d_val = conjectured_gamma_k_n(k=fixed_k, n=n, d=d)
    s_d_val = len(compute_Sn(S=S, n=n, group_op=op_wreath_prod))
    print("for n = {}, Z_d was {} and S_d was {}, ratio S_d/Z_d is {}".format(n, z_d_val, s_d_val, s_d_val/z_d_val))

#compare_func_of_Z_d_with_S_d(3, 8, 13, 3, 3, lambda x:x)
def compare_func_of_Z_d_with_S_d(d_upper, n_lower, n_upper, k_lower, k_upper, func):
    for k in range (k_lower, k_upper + 1):
      for n in range(n_lower, n_upper + 1):
        for d in range(3, d_upper + 1):
          z_d_val =  conjectured_gamma_k_n(k=k, n=n, d=d)
          s_d_val =  compute_approx_gamma_k_n(k = k, n = n, d = d, num_iters=20)
          print(f"n = {n}, k = {k}, d = {d} we get log(s_d/z_d) = {log(s_d_val/z_d_val)}")


# try num_iters different sets S with cardinality k, 
# compute |S^n| for each and take the max
# ex: compute_approx_gamma_k_n(k = 3, n = 8, d = 3, num_iters = 20, noisy = True)
def compute_approx_gamma_k_n(k, n, d, num_iters, noisy = False):
  best = -1
  best_set = set()

  s_d = list(make_s_d(d))

  def get_set(k, d):  
    S = set()
    # add something to the set k times 
    for _ in range(0, k):
      rand_s_d_elt = s_d[randint(0, factorial(d) - 1)]
      S.add((get_sparse_d_tuple(d), rand_s_d_elt))
    return S

  for _ in range(0, num_iters):
    S = get_set(k = k, d = d)
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

def test_S_d_bound(d_upper, n_lower, n_upper, k_lower, k_upper, bound_n_k_d):
    for k in range (k_lower, k_upper + 1):
      for n in range(n_lower, n_upper + 1):
        for d in range(3, d_upper + 1):
          bound_val =  bound_n_k_d(n, k, d)
          s_d_val =  compute_approx_gamma_k_n(k = k, n = n, d = d, num_iters=20)
          print(f"n = {n}, k = {k}, d = {d} we get s_d/bound = {s_d_val/bound_val} and s_d/(k^n) = {s_d_val/(k**n)}")

def test_S_d_bound_with_set(d_upper, n_lower, n_upper, k_lower, k_upper, bound_n_k_d, good_set):
    for k in range (k_lower, k_upper + 1):
      for n in range(n_lower, n_upper + 1):
        for d in range(3, d_upper + 1):
          bound_val =  bound_n_k_d(n, k, d)
          s_n = compute_Sn(good_set, n, op_wreath_prod)
          s_d_val =  len(s_n)
          print(f"n = {n}, k = {k}, d = {d} we get s_d/bound = {s_d_val/bound_val} and s_d/(k^n) = {s_d_val/(k**n)}")


#find_a_good_set(d = 3, k = 5, n = 6, num_iters=100)
def find_a_good_set(d, k, n, num_iters):
  best_set = set()
  best_val = 0
  prs = []

  #value_set = set()

  for _ in range(0, num_iters):
    S = get_rand_sparse_set(k=k, d=d)

    #num_distinct_second_components = len({elt[1] for elt in S})
    cur_val = len(compute_Sn(S = S, n = n, group_op= op_wreath_prod)) 
    
    #value_set.add(cur_val)
    #prs.append((cur_val, num_distinct_second_components))

    #print("With {} 2nd components we got {}".format(num_distinct_second_components, cur_val))
    if cur_val > best_val:
      best_set = S 
      best_val = cur_val 

  #print("we got {} different values".format(len(value_set)))
  print(best_val)
  #print(sorted(prs))
  #print_second_components(best_set)
  print(best_set)

#find_a_good_set(d = 3, k = 4, n = 10, num_iters=100)
good_set_ = {((58760462, 671770735, 368002213), (2, 1, 0)), ((831745905, 848100975, 907675504), (1, 0, 2)), ((821077691, 59173342, 933562559), (0, 1, 2)), ((660488124, 588224254, 236702212), (0, 2, 1))}

#something(5, 10, 2, 6, 3, 15)
def something(n_l, n_h, k_l, k_h, d, num_iters):
  for k in range (k_l, k_h + 1):
    for n in range(n_l, n_h + 1):
      s_d = compute_approx_gamma_k_n(k=k, n=n, d=d, num_iters=num_iters)
      z_d = conjectured_gamma_k_n(k=k, n=n, d=d)

      print("n = {}, k = {}, s_d/n^(d(k-1)) = {}".format(n, k, z_d/(n**(d*(k-1)))))

#something(5, 10, 2, 6, 3, 15)
# check if a set S is pairwise non_commutative
def is_pairwise_non_commutative(S):
  perms = [elt[1] for elt in S]
  k = len(perms)
  for i in range(0, k):
    for j in range(0, k):
      if i != j and op_s_d(perms[i], perms[j]) == op_s_d(perms[j], perms[i]):
        return False
  return True

#test_if_pairwise_noncommutativity_helps(6, 5, 3, 100)
def test_if_pairwise_noncommutativity_helps(n, k, d, num_iters):
  best_sets = []
  best_val = 0

  for _ in range(0, num_iters):
    S = get_rand_sparse_set(k=k, d=d)
    cur_val = len(compute_Sn(S=S, n=n, group_op=op_wreath_prod))

    if cur_val > best_val:
      best_sets = [S]
      best_val = cur_val
    elif cur_val == best_val:
      best_sets.append(S)
  
  non_com = 0
  com = 0
  for S in best_sets:
    print_second_components(S)
    if is_pairwise_non_commutative(S):
      non_com += 1
    else:
      com += 1
  
  print("Of the best sets, {} were non com and {} were com".format(non_com, com))


# return a set of tuples, where the size of this set is the same 
# as the size of S_n. In each tuple, the first compute is an element of s_n
# and the 2nd component is a list of lists, where each sublist is a list 
# of tags of the elements used to make that product

# actually, a map from the elements of s_n to this list of lists
def compute_s_n_and_how_we_achieved_each_elt(tagged_Set, n):
  if n == 1: 
    return {elt:[[tag]] for (tag, elt) in tagged_Set}
  
  one_less = compute_s_n_and_how_we_achieved_each_elt(tagged_Set, n - 1)

  return_dict = dict()
  for (tag, elt) in tagged_Set:
    for s_n_minus_one_elt in one_less.keys():
      potential_new_elt = op_wreath_prod(elt, s_n_minus_one_elt) 
      all_ways_to_get_here = []
      for way_to_get_s_n_minus_one_elt in one_less[s_n_minus_one_elt]:
        all_ways_to_get_here.append([tag] + way_to_get_s_n_minus_one_elt)

      if potential_new_elt not in return_dict:
        return_dict[potential_new_elt] = []
      return_dict[potential_new_elt] += all_ways_to_get_here
  return return_dict

def make_tagged_set(k, d):
  S = get_rand_sparse_set(k, d)
  tagged_Set = set()
  tag = 0
  for elt in S:
    tagged_Set.add((tag, elt))
    tag += 1
  return tagged_Set


#investigate_if_signature_determines_everything(11, 3, 3, lambda elt: (elt[0][0], elt[0][1]))
# signature is a function from Z_d prod S_d to anything
def investigate_if_signature_determines_everything(n, k, d, signature):
  S = get_rand_sparse_set(k, d)
  s_n = compute_Sn(S, n, op_wreath_prod)
  sig_set = {signature(elt) for elt in s_n}
  print(f"for n = {n}, k = {k}, d = {d} in this case, given signature", end = " ")
  if len(sig_set) == len(s_n):
    print("DID", end = " ")
  else:
    print("DID NOT", end = " ")
  print("uniquely determine the elements")
  print(f"|elts|/|sig| = {len(s_n) / len(sig_set)}")



def investigate_when_elts_same(n, k, d):
  def determine_common_elts(list_of_lists):
    common_len = len(list_of_lists[0])
    common_elt_list = []
    dif_elt_list = []
    for i in range(0, common_len):
      set_of_ith_elts = {lst[i] for lst in list_of_lists}
      if len(set_of_ith_elts) == 1:
        common_elt_list.append(i)
      else:
        dif_elt_list.append(i)
    return (common_elt_list, dif_elt_list)


  tagged_Set = make_tagged_set(k, d)
  second_comps = {tag: elt[1] for (tag, elt) in tagged_Set}
  elt_to_LL = compute_s_n_and_how_we_achieved_each_elt(tagged_Set, n)

  ways_to_get_to_element_to_frequency = {}


  for elt in elt_to_LL:
    ways_to_get_there = elt_to_LL[elt]

    # num_ways = len(ways_to_get_there)
    # if num_ways not in ways_to_get_to_element_to_frequency:
    #   ways_to_get_to_element_to_frequency[num_ways] = 0
    # ways_to_get_to_element_to_frequency[num_ways] += 1

    if len(ways_to_get_there) > 1:
      print(f"There are {len(ways_to_get_there)} ways to get this element:")
      for way in ways_to_get_there:
        print(way)
      #print(f"(shared, different) is {determine_common_elts(ways_to_get_there)}")
      print()

  print(f"n = {n}, k = {k}, d = {d}")
  print(ways_to_get_to_element_to_frequency)
  print(f'tag: second components are: {second_comps}')
  print(f'|S^n| = {len(elt_to_LL)}, and k^n = {k**n}')
  print(f'Just FYI, for Z^d prod Z_d the value is {conjectured_gamma_k_n(k = k, n = n, d = d)}')

#investigate_when_elts_same(10, 3, 3)
#investigate_when_elts_same(4, 3, 3)

def see_how_many_elts_from_multi_set(n, k, d, noisy = False):
  #s_d = list(make_s_d(d))
  #our_set = [(get_sparse_d_tuple(d), s_d[i]) for i in range(0, k)]
  our_set = list(get_rand_sparse_set(k, d))
  id_wreath = (tuple([0] * d), tuple([i for i in range(0, d)]))
  
  if noisy:
    print(f"Our ordered base set of elements is {our_set}")


  multisets = [list(p) for p in partition(n, k)]

  if len(multisets) != comb(n + k - 1, k - 1):
    raise ValueError("something went wrong, we don't have the right number of multisets")

  biggest = -1

  cur_sum = 0
  for ms in multisets:
    if noisy:
      print("Start of multiset")
      print(f"The multiset is {ms}")
    cur_set = make_all_products_from_multiset(ms, our_set, id_wreath, op_wreath_prod)
    ms_prod_size = len(cur_set)
    biggest = max(biggest,ms_prod_size)
    cur_sum += ms_prod_size

    if noisy:
      print(ms_prod_size)
  
  eventual_bound = comb(n + k - 1, k - 1) * biggest
  size_s_n = len(compute_Sn(our_set, n = n, group_op=op_wreath_prod))

  if size_s_n != cur_sum:
    print("Sanity check of |S^n| = sum over sets FAILED")
  else: 
    print("Sanity check of |S^n| = sum over sets GOOD")


  print_vals(n, k, d)
  print(f"We have |S^n| = {size_s_n}")
  print(f"Sum over the disjoint sets is {cur_sum}")
  print(f"To know how bad our bound will be, |S^n|/bound = {size_s_n/eventual_bound}")
  print(f"Largest number of products coming from one multiset was {biggest}")
  print(f"We have (d!)^n = {factorial(d)**(n-1)}")



#see_how_many_elts_from_multi_set(6, 6, 3)
#v = compute_approx_gamma_k_n(3, 10, 3, 200, True)
#print(comb(12, 2) * 2877)

  
#compute_s_n_and_what_comp_extracted_into_fst_snd(make_tagged_set(5, 3), 7)
def compute_s_n_and_what_comp_extracted_into_fst_snd(tagged_Set, n):
  tagsToElts = {tag:elt[1] for (tag, elt) in tagged_Set}
  # these lists should have length n
  def compute_jth_comp_list(tag_list, j):
    fst_comp_list = [0]
    cur_action = tagsToElts[tag_list[0]]
    for i in range(1, len(tag_list)):
      t = tag_list[i]
      fst_comp_list.append(cur_action[j])
      cur_action = op_s_d(cur_action, tagsToElts[t])
    return fst_comp_list

  elt_to_LL = compute_s_n_and_how_we_achieved_each_elt(tagged_Set, n)
  for elt in elt_to_LL:
    ways_to_get_here = elt_to_LL[elt]
    if len(ways_to_get_here) > 3:
      print("Start element:")
      for way in ways_to_get_here:
        fst_comps = compute_jth_comp_list(way, 0)
        snd_comps = compute_jth_comp_list(way, 1)
        print(f"first components: {fst_comps}, second components: {snd_comps}, tags: {way}")
      print()
  print(tagsToElts)

def compute_s_n_and_running_prods(tagged_Set, n, d):
  tagsToElts = {tag:elt[1] for (tag, elt) in tagged_Set}
  # these lists should have length n
  s_d = list(make_s_d(d))
  s_d_reverse_map = {s_d[i]:i for i in range(0, len(s_d))}

  def compute_running_s_d_products(tag_list):
    res = []
    cur_prod = tagsToElts[tag_list[0]]
    for i in range(0, len(tag_list)):
      t = tag_list[i]
      res.append(s_d_reverse_map[cur_prod])
      cur_prod = op_s_d(cur_prod, tagsToElts[t])
    return res

  elt_to_LL = compute_s_n_and_how_we_achieved_each_elt(tagged_Set, n)
  for elt in elt_to_LL:
    ways_to_get_here = elt_to_LL[elt]
    if len(ways_to_get_here) > 10:
      print("Start element:")
      # lets make a dictionary running products -> multiplicities
      prodToMul = {}
      for way in ways_to_get_here:
        running_prods = tuple(compute_running_s_d_products(way))
        if running_prods not in prodToMul:
          prodToMul[running_prods] = 0
        prodToMul[running_prods] += 1
      print(f"running product: {prodToMul}, tags: {way}")
      print()
  print(tagsToElts)

#do_stuff(8, 5, 3)
def do_stuff(n, k, d):
  compute_s_n_and_running_prods(make_tagged_set(k = k, d = d), n = n, d = d)

# we can't confidently test the exponent. In the Z_d case we only matched the 
# true exponent for k = 3 with n = 1000
def test_exponent(n, k, d, conj_exp):
  pass

# confirmed not a bound: n = 13, k = 4, d = 3
def bound(n, k, d):
  return comb(n + (k * d) - 1, (k * d) - 1)

# bound^2 is wayyyyy too big

# this is decent? must not be tight for all values of d though
def bound2(n, k, d):
  return abelian(n, k, d) * factorial(d) * (n-1)

def test_bound_2():
  for n in range (10, 100):
    for k in range (10, 100):
      for d in range(3, 10):
        print(f"n = {n}, d = {d}, k = {k} ratio k^n/bound2 is {k**n / bound2(n, k, d)}")

# from the paper: H < G then gamma_H(k, n) \le [G:H]^3 \gamma_G([G:H]^5 k, n)
def bound3(n, k, d):
  index = factorial(d - 1)
  C_1 = index**3
  C_2 = index**5
  return C_1 *conjectured_gamma_k_n(k= C_2 * k, n = n, d = d)


# try tighter constants than bound3
def bound4(n, k, d):
  index = factorial(d - 1)
  C_1 = index
  C_2 = index
  return conjectured_gamma_k_n(k= C_2 * k, n = n, d = d)

def bound5(n, k, d):
  return abelian(n, k, d) * factorial(d)**2 * n


#investigate_when_elts_same(10, 3, 3)
#test_S_d_bound(3, 9, 15, 3, 4, bound5)
#test_S_d_bound_with_set(3, 9, 14, 4, 4, bound2, good_set_)
#test_S_d_bound_with_set(3, 9, 13, 4, 4, bound, get_rand_sparse_set(4, 3))

# this is d = 3, k = 3
def get_transposition_set():
  one_with_two = (1, 0, 2)
  one_with_three = (2, 1, 0)
  two_with_three = (0, 2, 1)

  return {(get_sparse_d_tuple(3), one_with_two), (get_sparse_d_tuple(3), one_with_three), (get_sparse_d_tuple(3), two_with_three)}

# d = 4, k = 4
def get_klein_four_group_set():
  id = (0, 1, 2, 3)
  fst = (1, 0, 3, 2)
  snd = (2, 3, 0, 1)
  thd = (3, 2, 1, 0)

  return {(get_sparse_d_tuple(4), id), (get_sparse_d_tuple(4), fst), (get_sparse_d_tuple(4), snd), (get_sparse_d_tuple(4), thd)}


def test_klein_thing():
  #compute_approx_gamma_k_n(d = 4, k = 4, n = 10, num_iters= 15, noisy=True)
  print(len(compute_Sn(S = get_klein_four_group_set(), n = 12, group_op=op_wreath_prod)))

def get_transposition_set(k, d):
  switches = list(generate_all_transpositions_in_s_d(d))
  S = {(get_sparse_d_tuple(d), switches[i % comb(d, 2)]) for i in range(0, k)}
  return S


def see_if_transpositions_are_good_many_vals():
  def see_if_transpositions_are_good(n, k, d):
    iters = 20
    S = get_transposition_set(k, d)
    many_sets = compute_approx_gamma_k_n(k, n, d, iters)
    transposition_set_val = len(compute_Sn(S, n, op_wreath_prod))

    print_vals(n, k, d)
    # print(f'Trying {iters} different sets best is {many_sets}')
    # print(f'With transpositions: {transposition_set_val}')
    print(f"Transpositions are {transposition_set_val * 100 /many_sets}% of the best")
    print()

  loop(5, 7, 3, 5, 5, 8, see_if_transpositions_are_good)

def loop(_n, N, _k, K, _d, D, func):
  for n in range(_n, N + 1):
    for k in range(_k, K + 1):
      for d in range(_d, D + 1):
        func(n, k, d)

#see_if_transpositions_are_good_many_vals()

def compute_S_n_transpositions(n, k, d):
  return len(compute_Sn(get_transposition_set(k, d), n, op_wreath_prod))

def test_if_lots_of_identity_is_good(n, k, d, num_ids):
  S = {get_id_wreath(d) for _ in range(0, num_ids)}
  non_ids = get_rand_sparse_set(k-num_ids, d)
  S.update(non_ids)

  val_with_lots_ids = len(compute_Sn(S, n, op_wreath_prod))
  val_with_rand_sets = compute_approx_gamma_k_n(k, n, d, 50, True)

  print(f"ids: {val_with_lots_ids}")
  print(f"non ids: {val_with_rand_sets}")

#test_if_lots_of_identity_is_good(8, 4, 4, 2)

def test_transposition_set_bound():
  def bound_transposition_set(n, k, d):
    bound = d**n * abelian(n, k)
    actual = compute_S_n_transpositions(n, k, d)
    print_vals(n, k, d)
    print(bound/actual)

  for n in range(5, 10):
    for k in range(3, 5):
      for d in range(2, 5):
        bound_transposition_set(n,k,d)

#test_transposition_set_bound()

#test_klein_thing()


#print(len(compute_Sn(S = get_transposition_set(), n = 12, group_op=op_wreath_prod)))
#compute_approx_gamma_k_n(3, 15, 3, 10, True)

#compare_func_of_Z_d_with_S_d(3, 8, 15, 3, 3, lambda x: x)
