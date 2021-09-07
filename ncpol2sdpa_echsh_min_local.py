"""
In this script we compute lower bounds on H_min(AB|X=0,Y=0,E) for
devices constrained by some eCHSH_2 winning probability
"""


def ent(SDP):
    # Returns the entropy of the computed solution
    return -1 * log2(-SDP.dual)


from itertools import product
from math import log2, sqrt
import fractions

import ncpol2sdpa as ncp

k = 2
# for k = 2 

# VazVid Game k = 2
# X = Y = {0, 1}
# A = {0, 1, 2, 3}
# B = {0, 1, 2, 3}
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


# Global level of NPA relaxation
LEVEL = 2
# Maximum VazVid score for k = 2
WMAX = 0.5 + sqrt(2)/4

results = {}
WVazVids = [
    0.75, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, WMAX
]
for WVazVid in WVazVids:
    results[str(WVazVid)] = []

    A_config = [2**k, 2**k]
    B_config = [2**k, 2**k]
    # Measurement operators
    A = [Ax for Ax in ncp.generate_measurements(A_config, 'A')]
    B = [By for By in ncp.generate_measurements(B_config, 'B')]
    W = ncp.generate_operators('W', 2**k, hermitian=True)

    substitutions = {}
    moment_ineqs = []
    moment_eqs = []
    operator_eqs = []
    operator_ineqs = []
    localizing_monos = []  # op_eqs are processed last so need to add three Nones to end

    # Adding the constraints for the measurement operators
    substitutions.update(ncp.projective_measurement_constraints(A, B))

    # Defining the vazvid inequality for k = 2
    vazvid_expr = 0
    for (a, b, x, y) in product(range(2**k), range(2**k), [0, 1], [0, 1]):
        if a == (2**k-1) and b == (2**k-1):
            vazvid_expr += game_func(a, b, x, y) * (1 - sum(a for a in A[x])) \
                           * (1 - sum(b for b in B[y]))
        elif a == (2**k-1) and b != (2**k-1):
            vazvid_expr += game_func(a, b, x, y) * (1 - sum(a for a in A[x])) * B[y][b]
        elif a != (2**k-1) and b == (2**k-1):
            vazvid_expr += game_func(a, b, x, y) * A[x][a] * (1 - sum(b for b in B[y]))
        else:
            vazvid_expr += game_func(a, b, x, y) * A[x][a] * B[y][b]
    # divide by probability for each input
    vazvid_expr /= 4.0

    # constraint on the cglmp score
    w_exp = WVazVid
    score_con = [vazvid_expr - w_exp]

    # Commutation constraints for W (I think the projective msmt thing already includes ones for A,B)
    for w in W:
        for Ax in A:
            for a in Ax:
                substitutions.update({w * a: a * w})
        for By in B:
            for b in By:
                substitutions.update({w * b: b * w})

    # \sum W_a <= I_{R'}
    operator_ineqs += [1 - (sum(w for w in W))]
    # positivity constraints for W_a
    operator_ineqs += [w for w in W]
    # We must specify localizing mmonomials for the constraints of the
    # problem but by specifying None ncpol2sdpa uses a default set
    localizing_monos += [None] * (len(W) + 2)

    moment_equalities = moment_eqs[:]
    moment_inequalities = moment_ineqs[:] + score_con[:]
    operator_equalities = operator_eqs[:]
    operator_inequalities = operator_ineqs[:]

    # We now specify some extra monomials to include in the relaxation
    extra_monos = []
    for w in W:
        for Ax in A:
            for a in Ax:
                for By in B:
                    for b in By:
                        extra_monos += [a * b * w]
                extra_monos += [a * w]
        for By in B:
            for b in By:
                extra_monos += [b * w]

    result_sum = 0
    for x in range(2):
        # Objective function
        obj = sum(A[x][i] * W[i] for i in range(len(W)-1)) + \
            (1 - sum(a for a in A[x])) * W[len(W)-1]

        ops = ncp.flatten([A, B, W])
        sdp = ncp.SdpRelaxation(ops, verbose=1, normalized=True, parallel=0)
        sdp.get_relaxation(level=LEVEL,
                           equalities=operator_equalities,
                           inequalities=operator_inequalities,
                           momentequalities=moment_equalities,
                           momentinequalities=moment_inequalities,
                           objective=-obj,
                           substitutions=substitutions,
                           extramonomials=extra_monos,
                           localizing_monomials=localizing_monos)
        sdp.solve('mosek')
        print(
            f"For a vazid score {w_exp} and input x={x} we find a dual of {sdp.dual} "
            f"(primal {sdp.primal}) and with that an entropy of {ent(sdp)}.")
        result_sum += ent(sdp)
    # should divide result by 2.0, but we do that later when rendering
    results[str(WVazVid)] += [result_sum / 1.0]
print(results)
