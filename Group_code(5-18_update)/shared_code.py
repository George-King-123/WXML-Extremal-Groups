from random import randint
from math import comb

BIG_NUM = 10000000

def get_sparse_d_tuple(d):
    return tuple([randint(0, BIG_NUM) for _ in range(0, d)])

def compute_Sn(S, n, group_op):
  if n < 1:
    raise ValueError("should not give compute_Sn any value < 1")
  if n == 1:
    return S
  
  prod_with_one_less = compute_Sn(S, n - 1, group_op)
  s_n = set()
  for elt in S:
    for prod in prod_with_one_less:
      s_n.add(group_op(prod, elt))

  
  return s_n

# return a list of all tuples (x_1, ..., x_d)
# where the x_i are integers in {0, ..., k} such that
# sum x_i = k
def partition(k, d):
    if d == 1:
        return [(k,)]
    if k == 0:
        return [(0,) * d]
    partitions = []
    for i in range(k + 1):
        for sub_partition in partition(k - i, d - 1):
            partitions.append((i,) + sub_partition)
    return partitions

# set cur_prod to be the group identity when first called 
# multi_set is a list of length k s.t. multi_set[i] is the number 
# of occurences of S_list[i] in our multiset
# prods should be empty on first call
def make_all_products_from_multiset(multi_set, S_list, group_identity, group_op):
  k = len(S_list)

  #assert len(multi_set) == k

  all_prods = set()

  def helper(ms, cur_prod):
    if sum(ms) == 0:
      all_prods.add(cur_prod)
      return 
  
    for i in range(0, k):
      if ms[i] > 0:
        ms[i] -= 1
        new_prod = group_op(cur_prod, S_list[i])
        helper(ms, new_prod)
        ms[i] += 1

  helper(multi_set, group_identity)
  return all_prods


def abelian(n, k):
  return comb(n + k - 1, k - 1)

def print_vals(n, k, d):
   print(f"n = {n}, k = {k}, d = {d}")