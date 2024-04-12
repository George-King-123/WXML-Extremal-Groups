# TL;DR: Change line 62 to different values and run it

# in this file everything is 0 indexed.
# so f: [n] -> [k] is really f:{0, 1, ..., n-1} -> {0, 1, ..., k-1}

# the runtime of this simulation is (much) better than O(k^n), it is closer
# to the actual number of products. It has the same runtime as choosing sparse
# integers and computing products in recursive layers. 

# I can describe why it works if desired, but it's a little involved. 
# I did check that this gives the same results as the slow simulation, which clearly does work
# You can also check that it matches the formula for the abelian case when there are no negative signs

def compute_prods(S, n, k):
  # nested inner method
  # given a 
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
  len_n_minus_one_prods = compute_prods(S, n-1, k)
  len_n_prods = set()
  for prod in len_n_minus_one_prods:
    for i in range(0, k):
      len_n_prods.add(compute_one_more(prod, i, S))
  # last fcn_value is f(n) for our current value of n

  return len_n_prods

# requires num_neg <= k
def compute_size(n, k, num_neg):
  S = [1 for i in range (num_neg, k)] + [-1 for i in range(0, num_neg)]
  return len(compute_prods(S, n, k))

def run_sim(n, k):
  print("n = {} and k = {}".format(n, k))
  max_size = 0
  num_signs_max = 0
  for i in range(0, k+1):
    cardinality_s_n = compute_size(n, k, i)
    if cardinality_s_n > max_size:
      max_size = cardinality_s_n
      num_signs_max = i
    # Uncomment this to get actual sizes
    # print("|S^n| = {} when we have {} negative signs".format(cardinality_s_n, i))
  print("We achieved max size at {} signs".format(num_signs_max))


# change this to mess around with stuff
# when n > k it seems like fewer signs is better (like 1)
# when n = k it seems like ~2 is the best
# when n < k more seems better to have more signs

# you can make n, k larger than this, just let it run
# you can add print statements to the run_sim loop to check how much progress 
#   is being made
run_sim(n = 7, k = 7)
