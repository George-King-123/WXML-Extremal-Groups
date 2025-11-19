from random import randint

def incl_range(start, end):
  return range(start, end+1)

def get_sparse_d_tuple(d):
    BIG_NUM = 1000000000
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
      s_n.add(group_op(elt, prod))
  
  return s_n