"""
In this script we calculate the maximum CGLMP_3 win probability
"""

import ncpol2sdpa as ncp
from itertools import product
import fractions

# CGLMP game dimension 3
# X = Y = {0, 1}
# A = {0, 1, 2}
# B = {0, 1, 2}
# Probabilistic winning function
#  w(a,b,x,y,) = [a-b = xy mod 3] + 1/2 * [a-b = (-1)^{x \oplus y} + xy mod 3]
def game_pred(a, b, x, y):
    result = 0
    result += (a - b) % 3 == (x * y) % 3
    result += fractions.Fraction(((a - b) % 3 == (pow(-1, x ^ y) + x * y) % 3), 2)
    return result

level = 2
P = ncp.Probability([3, 3], [3, 3])
# CGLMP_3 win probability
CGLMP_game = 0
for (a, b, x, y) in product([0, 1, 2], [0, 1, 2], [0, 1], [0, 1]):
    CGLMP_game += game_pred(a, b, x, y) * P([a, b], [x, y])
objective = -CGLMP_game / 4
sdp = ncp.SdpRelaxation(P.get_all_operators())
sdp.get_relaxation(level, objective=objective,
                   substitutions=P.substitutions)
sdp.solve('mosek')
print(-sdp.primal)
print(-sdp.dual)
