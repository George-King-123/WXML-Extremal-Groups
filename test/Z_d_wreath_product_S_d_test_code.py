from shared_code import get_sparse_d_tuple, compute_Sn
from group_operations import op_wreath_prod
import Z_d_product_Z_d

# Since Z_d is a special case of S_d, we can test the S_d code by making sure a set 
# that only has "shift by one" permutations agrees with the formula we came up 
# with for Z_d
# using the given n, k, d, and t many components with 1 in the 2nd component
# with k - t elements having 0 in the 2nd component
def restrict_to_Z_d(n, k, d, t):
  def get_s_d_set():
    # shift_by_one \in Z_d
    shift_by_one = tuple([i for i in range(1, d)] + [0])
    identity = tuple([i for i in range(0, d)])
    
    S = set()
    
    # add something to the set k times 
    for i in range(0, k):
      if i < t:
        S.add((get_sparse_d_tuple(d), shift_by_one))
      else:
        S.add((get_sparse_d_tuple(d), identity))
    return S
  
  val_from_s_d_sim = len(compute_Sn(S = get_s_d_set(), n = n, group_op = op_wreath_prod))
  val_from_z_d_formula = Z_d_product_Z_d.formula_when_all_zero_one_in_snd_comp(n = n, k =k, p=t, d=d)

  return val_from_s_d_sim == val_from_z_d_formula

# ex: test_a_bunch_of_restrictions(6, 6, 3)
def test_a_bunch_of_restrictions(n_upper, k_upper, d_upper):
  for n in range(1, n_upper + 1):
    for k in range(1, k_upper + 1):
      for d in range(1, d_upper + 1):
        for t in range(0, k + 1):
          if not restrict_to_Z_d(n=n, k = k, d= d, t = t):
            print("fail")

  print("success")