from random import randint
from math import comb as choose_built_in

BIG_NUM = 1_000_000_000

def big_rand_num():
  return randint(0, BIG_NUM)

def get_sparse_d_tuple(d):
    return tuple([big_rand_num() for _ in range(0, d)])

def choose(n, k):
  if (k < 0 or n < 0):
    return 0
  return choose_built_in(n, k)

def format_large_num(num): 
  return f"{num:,}"

def incl_range(start, end):
  return range(start, end+1)

def compute_Sn(S, n, group_op):
  if n < 1:
    raise ValueError("should not give compute_Sn any value < 1")
  if n == 1:
    return S
  
  prod_with_one_less = compute_Sn(S, n - 1, group_op)
  s_n = set()
  for elt in S:
    for prod in prod_with_one_less:
      s_n.add(group_op(elt, prod))
  
  return s_n