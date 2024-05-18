import D_infinity, Z_d_product_Z_d, Z_d_wreath_product_S_d
from group_operations import op_Z_d_Z_d, op_d_inf, op_wreath_prod
from shared_code import compute_Sn, make_all_products_from_multiset, partition, abelian, BIG_NUM
from random import randint

# PRINT_EVERY_ASSERTION = True

# def my_assert(cond):
#   if PRINT_EVERY_ASSERTION:
#     if not cond:
#       print("fail")
#     else:
#       print("success")
#   assert not cond

# test make_all_products_from_multiset
def test_multiset():
  def get_multi_sets(n, k):
    multi_sets = [list(p) for p in partition(n, k)]
    assert len(multi_sets) == abelian(n, k)
    return multi_sets


  # test on the integers, why not. We should only ever get 1
  def test_integers(n, k):
    our_set = [(n+1) ** i for i in range(1, k+1)]
    multi_sets = get_multi_sets(n, k)

    for ms in multi_sets:
      all_prods = make_all_products_from_multiset(ms, our_set, 0, lambda x, y: x + y)
      assert len(all_prods) == 1
  
  #test_integers(5, 5)

  # no real way to test this on D_inf it turns out. We don't have that |S^n| is a 
  # sum over the multi-sets because of cancellation: we can end up with the same product
  # from two different multi-sets. 
  # def test_on_D_inf(n, k, t):
  #   our_set = [(randint(0, BIG_NUM),1) for i in range (t, k)] + [(randint(0, BIG_NUM),-1) for i in range(0, t)]
  #   multi_sets = get_multi_sets(n, k)
  #   total = 0
  #   for ms in multi_sets:
  #     all_prods = make_all_products_from_multiset(ms, our_set, (0, 1), op_d_inf)
  #     total += len(all_prods)

  #   assert total == D_infinity.compute_s_n_with_formula(n, k, t)

  def test_on_Z_d(n, k, t, d):
    def id_Z_d(d):
      return (tuple([0] * d), 0)
    # make a set with k-t elements have 0 in the second component, 
    # t elements having a 1 in the second component
    our_set = list(Z_d_product_Z_d.make_S(k, d, tuple([k-t, t] + [0] * (d-2))))
    multi_sets = get_multi_sets(n, k)
    total = 0
    for ms in multi_sets:
      all_prods = make_all_products_from_multiset(ms, our_set, id_Z_d(d), op_Z_d_Z_d)
      total += len(all_prods)

    assert total == Z_d_product_Z_d.formula_when_all_one(n, k, t, d)

  # this passes
  # test_Z_d_many_vals(5, 7, 3, 5, 2, 4) 

  def test_Z_d_many_vals(n_l, n_u, k_l, k_u, d_l, d_u):
    for n in range(n_l, n_u + 1):
      for k in range(k_l, k_u + 1):
        for d in range(d_l, d_u + 1):
          for t in range(0, k+1):
            test_on_Z_d(n, k, t, d)

  # FAILS n = 6, k = 5, d = 4
  # 14457 != 14466
  # FAILS n = 6, k = 4, d = 4
  def test_on_S_d(n, k, d):
    def id_S_d(d):
      return (tuple([0] * d), tuple([i for i in range(0, d)]))
    # make a set with k-t elements have 0 in the second component, 
    # t elements having a 1 in the second component
    S = Z_d_wreath_product_S_d.get_rand_sparse_set(k, d)
    # print(S)
    multi_sets = get_multi_sets(n, k)

    our_s_n = set()
    our_cardinality_s_n = 0
    for ms in multi_sets:
      all_prods = make_all_products_from_multiset(ms, list(S), id_S_d(d), op_wreath_prod)
      our_s_n.update(all_prods)
      our_cardinality_s_n += len(all_prods)

    # checks that the unions were disjoint 
    assert our_cardinality_s_n == len(our_s_n)

    real_s_n = compute_Sn(S, n, op_wreath_prod)
    #print(our_s_n)
    #print(real_s_n)
    # print(f"real - actual: {real_s_n.difference(our_s_n)}")
    # print(f"actual - real: {our_s_n.difference(real_s_n)}")
    # print(len(real_s_n.difference(our_s_n)))
    # print(len(our_s_n.difference(real_s_n)))
    # print(our_cardinality_s_n)
    # print(len(real_s_n))
    # print("\n")



    assert len(real_s_n.difference(our_s_n)) == 0
    assert len(our_s_n.difference(real_s_n)) == 0
    assert our_cardinality_s_n == len(real_s_n)

  # this passes
  # test_S_d_many_vals(12, 12, 4, 4, 3, 3) 
  # test_S_d_many_vals(9, 11, 4, 4, 3, 3) 
  def test_S_d_many_vals(n_l, n_u, k_l, k_u, d_l, d_u):
    for n in range(n_l, n_u + 1):
      for k in range(k_l, k_u + 1):
        for d in range(d_l, d_u + 1):
          test_on_S_d(n, k, d)
  test_S_d_many_vals(3, 8, 3, 5, 4, 7) 

  
test_multiset()
  

