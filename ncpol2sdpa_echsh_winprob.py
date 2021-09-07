import ncpol2sdpa as ncp
from itertools import product
import fractions

# Results:
#  For k=2 just recover the CHSH win probability 0.5 + sqrt(2)/4

# number of parallel CHSH games
k = 2

# VazVid_k
# X = Y = {0, 1}
# A = {0, 1, ..., 2^k-1}
# B = {0, 1, ..., 2^k-1}
def game_func(a, b, x, y):
    if x == 0 and y == 0:
        return 1 if a == b else 0
    # compute relative Hammond distance
    a_bits = "{0:b}".format(a).rjust(k, "0")
    b_bits = "{0:b}".format(b).rjust(k, "0")
    n_diff = 0
    for i in range(len(a_bits)):
        n_diff += a_bits[i] != b_bits[i]
    dist = fractions.Fraction(n_diff, k)
    if y == 1:  # (B, 0)
        return 1 if dist <= fractions.Fraction(16, 100) else 0
    elif x == 1 and y == 0:
        return 1 if dist >= fractions.Fraction(49, 100) and dist <= fractions.Fraction(51, 100) else 0
    raise ValueError("Invalid x={x} and/or y={y}")

level = 2
P = ncp.Probability([2**k, 2**k], [2**k, 2**k])
# VazVid_k win probability
VazVid_game = 0
for (a, b, x, y) in product(range(2**k), range(2**k), [0, 1], [0, 1]):
    VazVid_game += game_func(a, b, x, y) * P([a, b], [x, y])
objective = -VazVid_game / 4
sdp = ncp.SdpRelaxation(P.get_all_operators(), verbose=0)
sdp.get_relaxation(level, objective=objective,
                   substitutions=P.substitutions)
sdp.solve('mosek')
print(-sdp.primal)
print(-sdp.dual)
