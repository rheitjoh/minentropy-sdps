"""
In this script we calculate H_min(A|X) for the CHSH game constrained by some
CHSH winning probability.
"""

def ent(SDP):
    # Returns the entropy of the computed solution
    return -1 * log2(-SDP.dual)


from math import sqrt, log2

import ncpol2sdpa as ncp

# Global level of NPA relaxation
LEVEL = 2
# Maximum CHSH score
WMAX = 0.5 + sqrt(2) / 4

results = {}
# the different CHSH win probabilities we want to calculate local min-entropy for
WCHSHs = [
    0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.805, 0.81,
    0.815, 0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, WMAX
]
for WCHSH in WCHSHs:
    results[str(WCHSH)] = []
    A_config = [2, 2]
    B_config = [2, 2]
    # Measurement operators
    A = [Ax for Ax in ncp.generate_measurements(A_config, 'A')]
    B = [By for By in ncp.generate_measurements(B_config, 'B')]
    W = ncp.generate_operators('W', 2, hermitian=True)

    substitutions = {}
    moment_ineqs = []
    moment_eqs = []
    operator_eqs = []
    operator_ineqs = []
    localizing_monos = []  # op_eqs are processed last so need to add three Nones to end

    # Adding the constraints for the measurement operators
    substitutions.update(ncp.projective_measurement_constraints(A, B))

    # Defining the chsh inequality
    chsh_expr = (A[0][0] * B[0][0] + (1 - A[0][0]) * (1 - B[0][0]) + \
                 A[0][0] * B[1][0] + (1 - A[0][0]) * (1 - B[1][0]) + \
                 A[1][0] * B[0][0] + (1 - A[1][0]) * (1 - B[0][0]) + \
                 A[1][0] * (1 - B[1][0]) + (1 - A[1][0]) * B[1][0]) / 4.0

    # constraint on the chsh score
    w_exp = WCHSH
    score_con = [chsh_expr - w_exp]

    # Commutation constraints for W (I think the projective msmt thing already includes ones for A,B)
    for w in W:
        for Ax in A:
            for a in Ax:
                substitutions.update({w * a: a * w})
        for By in B:
            for b in By:
                substitutions.update({w * b: b * w})

    # \sum W_{a,b} <= I_{R'}
    operator_ineqs += [1 - (W[0] + W[1])]
    # positivity constraints for W_{a,b}
    operator_ineqs += [W[0], W[1]]
    # We must specify localizing mmonomials for the constraints of the
    # problem but by specifying None ncpol2sdpa uses a default set
    localizing_monos += [None, None, None, None]

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
        for y in range(2):
            # Objective function
            obj = A[x][0] * W[0] + \
                  (1 - A[x][0]) * W[1]

            ops = ncp.flatten([A, B, W])
            sdp = ncp.SdpRelaxation(ops, verbose=0, normalized=True, parallel=0)
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
                f"For a chsh score {w_exp} and inputs x {x} y {y} we find an sdp dual value of {sdp.dual} and with that an entropy of {ent(sdp)}.")
            result_sum += ent(sdp)
    results[str(WCHSH)] += [result_sum / 4.0]
print(results)
