from math import sqrt

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import entropy


def draw_chsh():
    chsh = {'0.75': [5.364178620785489e-09], '0.76': [0.03351186500978613],
            '0.77': [0.0732431789849653], '0.78': [0.11871745961231302],
            '0.79': [0.17091135494736448], '0.8': [0.23129953462909042],
            '0.805': [0.2652622107464345], '0.81': [0.30226806580240667],
            '0.815': [0.3428757955108997], '0.82': [0.3878425023544511],
            '0.825': [0.43823252399993545], '0.83': [0.49562245193907073],
            '0.835': [0.5625231349300338], '0.84': [0.6433537160372045],
            '0.845': [0.7474134095392008], '0.85': [0.902766883413926],
            '0.8535533905932737': [1.2280672521048372]}
    chsh_local = {'0.75': [4.736578431240713e-09], '0.76': [0.030374687513070642],
                  '0.77': [0.06415028547758023], '0.78': [0.10199960647159066],
                  '0.79': [0.14484749545721354], '0.8': [0.19402122528358917],
                  '0.805': [0.22157083859643106], '0.81': [0.2515387105898424],
                  '0.815': [0.2843843521656115], '0.82': [0.3207276198225399],
                  '0.825': [0.36143815333710216], '0.83': [0.40780073230425873],
                  '0.835': [0.461853032244506], '0.84': [0.5271852906084178],
                  '0.845': [0.6113226569866042], '0.85': [0.7369655033151004],
                  '0.8535533905932737': [0.9997243499206518]}
    chsh_analytic = {'0.75': [], '0.76': [], '0.77': [], '0.78': [], '0.79': [], '0.8': [],
                     '0.805': [], '0.81': [], '0.815': [], '0.82': [], '0.825': [], '0.83': [],
                     '0.835': [], '0.84': [], '0.845': [], '0.85': [], '0.8535533905932737': []}
    for k in chsh_analytic.keys():
        chsh_analytic[k] = 1 - entropy([
            1 / 2 + 1 / 2 * (sqrt(16 * float(k) * (float(k) - 1) + 3)),
            1 - (1 / 2 + 1 / 2 * (sqrt(16 * float(k) * (float(k) - 1) + 3)))
        ], base=2)

    fig, ax = plt.subplots()
    ax.plot([k for k in chsh.keys()], [v[0] for v in chsh.values()], "b.",
            label=r"$H_\mathrm{min}(AB|E)$")
    ax.plot([k for k in chsh_local.keys()], [v[0] for v in chsh_local.values()], "g.",
            label=r"$H_\mathrm{min}(A|E)$")
    ax.plot([k for k in chsh_analytic.keys()], [v for v in chsh_analytic.values()], "r.",
            label=r"$H(A|E)$ analytic")
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MaxNLocator(8))
    plt.xlabel("CHSH win probability")
    plt.ylabel("Entropy in bits")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    fig.savefig(
        "/home/raphael/Documents/papers/masterthesis/thesis/figures/chsh_entropy_results.pdf",
        bbox_inches="tight")


def draw_cglmp_3():
    cglmp_local = {'0.75': [-7.991858783275906e-10], '0.76': [0.03041865821532462],
                   '0.77': [0.0674501062337608], '0.78': [0.11056454281956948],
                   '0.79': [0.1602725246698503], '0.8': [0.2177397784295937],
                   '0.805': [0.2499493775661259], '0.81': [0.2848964203947455],
                   '0.815': [0.32299873401514967], '0.82': [0.36479429361331517],
                   '0.825': [0.4109931831373559], '0.83': [0.4625612865027591],
                   '0.835': [0.5208608306635041], '0.84': [0.5879114966044465],
                   '0.845': [0.6669169508008348], '0.85': [0.7634852748041835],
                   '0.86': [1.0769211981590643], '0.8643567588466105': [1.5837128088186747]}

    fig, ax = plt.subplots()
    ax.plot([k for k in cglmp_local.keys()], [v[0] for v in cglmp_local.values()], "r.",
            label=r"$H_\mathrm{min}(A|E)$")
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MaxNLocator(8))
    plt.xlabel("CGLMP win probability")
    plt.ylabel("Entropy in bits")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    fig.savefig(
        "/home/raphael/Documents/papers/masterthesis/thesis/figures/cglmp_3_entropy_results.pdf",
        bbox_inches="tight"
    )

def draw_vazvid_2():
    vazvid_local = dict((w, e[0]/2) for w, e in {
        '0.75': [3.088883184243907e-10], '0.8': [0.3880419793037799], 
        '0.81': [0.5030767117733151], '0.82': [0.6414543513465181], 
        '0.83': [0.8156003064239574], '0.84': [1.054367346289726], 
        '0.85': [1.4739298818965756], 
        '0.8535533905932737': [1.9971284232041318]
    }.items())

    fig, ax = plt.subplots()
    ax.plot([k for k in vazvid_local.keys()], [v for v in vazvid_local.values()], "r.",
            label=r"$H_\mathrm{min}(A|E)$")
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MaxNLocator(7))
    plt.xlabel(r"$\mathrm{eCHSH}_2$ win probability")
    plt.ylabel("Entropy in bits")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    fig.savefig(
        "/home/raphael/Documents/papers/masterthesis/thesis/figures/eCHSH_2_entropy_results.pdf",
        bbox_inches="tight"
    )


if __name__ == "__main__":
    draw_vazvid_2()
