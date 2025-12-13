# compute the group operation in Z^d product S_d given two ordered pairs g and h
# consider elements of S_d to be bijections on {0, 1, ..., d-1}
def op_wreath_prod (g, h):
  perm = g[1]
  d = len(g[0])
  res = [0] * d
  for i in range(0, d):
    res[i] = g[0][i] + h[0][perm[i]]
  return (tuple(res), op_s_d(g[1], h[1]))

# again consider elements of s_d to be bijections on {0, 1, ..., d-1}
# given two d-tuples representing such bijections, compute the product 
# and return as a tuple
def op_s_d(sigma, tau):
  d = len(sigma)
  res = [0] * d
  for i in range(0, d):
    res[i] = sigma[tau[i]]
  return tuple(res)

# compute the group operation in Z^d product Z_d given two ordered pairs g and h
# g[0] and h[0] should be tuples of length d (i.e. the same length)
# g[1] and h[1] should be integers in {0, 1, ..., d-1} 
def op_Z_d_Z_d (g, h):
  g_fst = g[0]
  rot = g[1]
  h_fst = h[0]
  d = len(g_fst)
  res = [0] * d
  for i in range(0, d):
    # yes, you need the +d, an artifact of how % works in python
    res[i] = g_fst[i] + h_fst[(i - rot + d) % d]
  return (tuple(res), (g[1] + h[1]) % d)

# compute the group operation in D_infinity given two ordered pairs g and h
# g[0] and h[0] should be sparse integers, g[1] and h[1] should be signs, +1 or -1
def op_d_inf (g, h):
  return (g[0] + g[1] * h[0], g[1] * h[1])
