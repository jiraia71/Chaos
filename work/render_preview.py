"""Preview das figuras KAM em matplotlib (o mestre continua sendo fig_qzs_kam.m).
Reproduz K7 (mapas FLI 2.5 GL), K8 (secoes de Poincare 2.5 GL) e K9
(amortecimento 2.5 GL) a partir de dados_kam.mat, para inspecao rapida sem MATLAB.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.io import loadmat

MAT = sys.argv[1] if len(sys.argv) > 1 else "/home/user/Chaos/dados_kam.mat"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/user/Chaos/work"

NAMES = ["monostable", "shallow_wells", "deep_wells"]
LABELS = ["Monostable  (η=η_QZS)", "Shallow wells  (η=0.60)", "Deep wells  (η=0.30)"]


def cmap_from(hexes):
    return LinearSegmentedColormap.from_list("c", ["#" + h for h in hexes])


FLI = cmap_from(["0b1d3a", "1f4e8c", "3a8fb7", "9ad0c2", "f4e285", "f4a259", "bc4b51", "5b1a18"])
ENE = cmap_from(["2b1a6f", "3b5bdb", "1c9fd6", "20b2aa", "6cc24a", "e0c300", "f08c00", "d7263d"])


def unwrap(x):
    while isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
        x = x.reshape(-1)[0]
    return x


def U(x, eta):
    return 1.55 * x * x - 3 * np.sqrt((1.5 * eta) ** 2 + x * x)


def load():
    D = unwrap(loadmat(MAT, squeeze_me=False, struct_as_record=False)["D"])
    return {nm: unwrap(getattr(D, nm)) for nm in NAMES}


def g(C, k):
    return np.asarray(getattr(C, k))


def gs(C, k):
    return float(np.asarray(getattr(C, k)).ravel()[0])


def separatrix(ax, C, XW, VW):
    xw = gs(C, "xw"); vw = gs(C, "vw"); eta = gs(C, "eta")
    Um = gs(C, "Umin"); depth = gs(C, "depth")
    xx = np.linspace(-xw, xw, 400); vv = np.linspace(-vw, vw, 400)
    XG, VG = np.meshgrid(xx, vv)
    HG = .5 * VG ** 2 + U(XG, eta) - Um
    if depth > 0:
        ax.contour(XG, VG, HG, [depth], colors="k", linewidths=.5, linestyles="--")


def fig_K7(D):
    fa = g(D["monostable"], "abs_f").ravel()
    fig, axes = plt.subplots(3, len(fa), figsize=(9, 8.2), constrained_layout=True)
    for r, nm in enumerate(NAMES):
        C = D[nm]; xs = g(C, "abs_xs").ravel(); vs = g(C, "abs_vs").ravel()
        eta = gs(C, "eta"); Um = gs(C, "Umin")
        maps = g(C, "abs_maps"); reg = g(C, "abs_regular").ravel()
        XG, VG = np.meshgrid(xs, vs); H0 = .5 * VG ** 2 + U(XG, eta) - Um
        for c in range(len(fa)):
            ax = axes[r, c]
            im = ax.imshow(np.clip(maps[c], 0, 20), origin="lower", aspect="auto",
                           extent=[xs[0], xs[-1], vs[0], vs[-1]], cmap=FLI, vmin=0, vmax=20)
            separatrix(ax, C, None, None)
            ax.contour(XG, VG, H0, [0.5], colors="w", linewidths=.4, linestyles=":")
            frac = float(np.mean(maps[c][H0 <= 0.5] <= 8))
            ax.text(.03, .96, f"regular {100*frac:.0f}%", transform=ax.transAxes,
                    va="top", color="w", fontsize=8)
            if r == 0:
                ax.set_title(f"f = {fa[c]:g}", fontsize=10)
            if c == 0:
                ax.set_ylabel(LABELS[r] + "\n dX/dt", fontsize=8)
            if r == 2:
                ax.set_xlabel("X₀", fontsize=9)
    fig.colorbar(im, ax=axes, shrink=.5, label="FLI após 400 períodos (escuro: toros KAM; claro: caos)")
    fig.suptitle("K7 — FLI do sistema completo QZS-ADV (2.5 GL), limite conservativo", fontsize=12)
    p = f"{OUT}/K7_FLI_2p5GL_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


def fig_K8(D):
    fs = g(D["monostable"], "sec2_f").ravel()
    fig, axes = plt.subplots(3, len(fs), figsize=(10, 8), constrained_layout=True)
    for r, nm in enumerate(NAMES):
        C = D[nm]; xw = gs(C, "xw"); vw = gs(C, "vw")
        pts = g(C, "sec2_pts"); fli = g(C, "sec2_fli"); E0 = g(C, "sec2_E0").ravel()
        for c in range(len(fs)):
            ax = axes[r, c]
            P = pts[c]  # (1500,2,64)
            chaotic = fli[c] > 10
            if chaotic.any():
                xg = P[:, 0, chaotic].ravel(); yg = P[:, 1, chaotic].ravel()
                ax.scatter(xg, yg, s=.2, c=[[.55, .55, .55]], linewidths=0)
            keep = ~chaotic
            if keep.any():
                xr = P[:, 0, keep]; yr = P[:, 1, keep]
                cc = np.clip(E0[keep] / .8, 0, 1)
                cc = np.repeat(cc[None, :], xr.shape[0], axis=0)
                ax.scatter(xr.ravel(), yr.ravel(), s=.3, c=cc.ravel(), cmap=ENE, vmin=0, vmax=1, linewidths=0)
            separatrix(ax, C, None, None)
            ax.set_xlim(-xw, xw); ax.set_ylim(-vw, vw)
            ax.text(.03, .96, f"caótico {100*np.mean(chaotic):.0f}%", transform=ax.transAxes,
                    va="top", fontsize=7, bbox=dict(fc="w", ec="none", pad=1))
            if r == 0:
                ax.set_title(f"f = {fs[c]:g}", fontsize=10)
            if c == 0:
                ax.set_ylabel(LABELS[r] + "\n dX/dt", fontsize=8)
            if r == 2:
                ax.set_xlabel("X", fontsize=9)
    fig.suptitle("K8 — Seções de Poincaré do QZS-ADV (2.5 GL); cinza: caótico (FLI>10); cor: energia inicial", fontsize=11)
    p = f"{OUT}/K8_secoes_2p5GL_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


def zlab(z1, z2):
    s = lambda z: "0" if z == 0 else f"1e{int(round(np.log10(z)))}"
    return f"ζ₁={s(z1)}, ζ₂={s(z2)}"


def fig_K9(D):
    z1 = g(D["monostable"], "damp2_zeta1").ravel(); z2 = g(D["monostable"], "damp2_zeta2").ravel()
    n = len(z1)
    fig, axes = plt.subplots(3, n, figsize=(10, 8), constrained_layout=True)
    for r, nm in enumerate(NAMES):
        C = D[nm]; xw = gs(C, "xw"); vw = gs(C, "vw")
        early = g(C, "damp2_early"); late = g(C, "damp2_late"); E0 = g(C, "damp2_E0").ravel()
        for c in range(n):
            ax = axes[r, c]
            e = early[c]; l = late[c]
            ax.scatter(e[:, 0, :].ravel(), e[:, 1, :].ravel(), s=.2, c=[[.81, .81, .81]], linewidths=0)
            cc = np.clip(E0 / .8, 0, 1); cc = np.repeat(cc[None, :], l.shape[0], axis=0)
            sz = .3 if (z1[c] == 0 and z2[c] == 0) else 3
            ax.scatter(l[:, 0, :].ravel(), l[:, 1, :].ravel(), s=sz, c=cc.ravel(), cmap=ENE, vmin=0, vmax=1, linewidths=0)
            separatrix(ax, C, None, None)
            ax.set_xlim(-xw, xw); ax.set_ylim(-vw, vw)
            if r == 0:
                ax.set_title(zlab(z1[c], z2[c]), fontsize=9)
            if c == 0:
                ax.set_ylabel(LABELS[r] + "\n dX/dt", fontsize=8)
            if r == 2:
                ax.set_xlabel("X", fontsize=9)
    fig.suptitle("K9 — Amortecimento no QZS-ADV (2.5 GL); cinza: 1os 150 per.; cor: per. 2701-3000 (atratores)", fontsize=11)
    p = f"{OUT}/K9_amortecimento_2p5GL_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


if __name__ == "__main__":
    D = load()
    for fn in (fig_K7, fig_K8, fig_K9):
        print(fn(D), flush=True)
