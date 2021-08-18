"""
In this script we calculate the maximum CHSH win probability
"""

import ncpol2sdpa as ncp

level = 2
P = ncp.Probability([2, 2], [2, 2])
# CHSH win probability
CHSH_game = P([0, 0], [0, 0]) + P([1, 1], [0, 0]) + \
            P([0, 0], [1, 0]) + P([1, 1], [1, 0]) + \
            P([0, 0], [0, 1]) + P([1, 1], [0, 1]) + \
            P([0, 1], [1, 1]) + P([1, 0], [1, 1])
objective = -CHSH_game/4
sdp = ncp.SdpRelaxation(P.get_all_operators())
sdp.get_relaxation(level, objective=objective,
                   substitutions=P.substitutions)
sdp.solve('mosek')
print(-sdp.primal)
print(-sdp.dual)
