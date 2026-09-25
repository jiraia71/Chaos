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


def qcolor(q):
    return {1: (.84, .15, .24), 2: (.94, .55, 0), 3: (.18, .62, .27)}.get(int(q), (.11, .49, .84))


def fig_K1(D):
    fig, axes = plt.subplots(2, 3, figsize=(11, 6), constrained_layout=True)
    for c, nm in enumerate(NAMES):
        C = D[nm]; eta = gs(C, "eta"); xw = gs(C, "xw"); Um = gs(C, "Umin"); depth = gs(C, "depth")
        # topo: potencial + ressonâncias
        ax = axes[0, c]
        x = np.linspace(-xw, xw, 1200); u = U(x, eta) - Um
        ax.plot(x, u, color=(.13, .13, .13), lw=1.1)
        if depth > 0:
            ax.axhline(depth, ls="--", color=(.4, .4, .4), lw=.6)
        rdE = g(C, "res_dE").ravel(); rq = g(C, "res_q").ravel()
        for e, q in zip(rdE, rq):
            xx = x.copy(); yy = np.full_like(xx, e); yy[u > e] = np.nan
            ax.plot(xx, yy, color=qcolor(q), lw=.9)
        ax.set_xlim(-xw, xw); ax.set_ylim(-.03, .95); ax.set_title(LABELS[c], fontsize=9)
        ax.set_xlabel("X")
        if c == 0:
            ax.set_ylabel("U(X) − U_min")
        # base: omega(E) por ramo
        ax = axes[1, c]
        dE = g(C, "dE").ravel(); w = g(C, "omega_quad").ravel(); over = g(C, "over").ravel()
        for b in (0, 1):
            m = over == b
            if m.sum() < 3:
                continue
            o = np.argsort(dE[m]); col = (.11, .25, .48) if b == 1 else (.70, .15, .12)
            ax.plot(dE[m][o], w[m][o], color=col, lw=1.2)
        rrho = g(C, "res_rho").ravel()
        for e, q, rho in zip(rdE, rq, rrho):
            ax.plot([e], [rho], "o", ms=4, mfc=qcolor(q), mec="w", mew=.4)
            ax.axhline(rho, ls=":", color=qcolor(q), lw=.4)
        ax.set_xscale("log"); ax.set_xlim(1e-5, 1.2); ax.set_ylim(0, 1.75)
        ax.set_xlabel("E − U_min")
        if c == 0:
            ax.set_ylabel("ω(E)/Ω")
    fig.suptitle("K1 — Esqueleto integrável (f=0): potencial, ressonâncias p/q e torção dω/dE", fontsize=12)
    p = f"{OUT}/K1_esqueleto_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


def _sections(D, ptsk, flik, e0k, fk, title, fname):
    fs = g(D["monostable"], fk).ravel()
    fig, axes = plt.subplots(3, len(fs), figsize=(10, 8), constrained_layout=True)
    for r, nm in enumerate(NAMES):
        C = D[nm]; xw = gs(C, "xw"); vw = gs(C, "vw")
        pts = g(C, ptsk); fli = g(C, flik); E0 = g(C, e0k).ravel()
        for c in range(len(fs)):
            ax = axes[r, c]; P = pts[c]; chaotic = fli[c] > 10
            if chaotic.any():
                ax.scatter(P[:, 0, chaotic].ravel(), P[:, 1, chaotic].ravel(), s=.2, c=[[.55, .55, .55]], linewidths=0)
            keep = ~chaotic
            if keep.any():
                xr = P[:, 0, keep]; yr = P[:, 1, keep]
                cc = np.repeat(np.clip(E0[keep] / .8, 0, 1)[None, :], xr.shape[0], axis=0)
                ax.scatter(xr.ravel(), yr.ravel(), s=.3, c=cc.ravel(), cmap=ENE, vmin=0, vmax=1, linewidths=0)
            separatrix(ax, C, None, None); ax.set_xlim(-xw, xw); ax.set_ylim(-vw, vw)
            ax.text(.03, .96, f"caótico {100*np.mean(chaotic):.0f}%", transform=ax.transAxes, va="top", fontsize=7,
                    bbox=dict(fc="w", ec="none", pad=1))
            if r == 0:
                ax.set_title(f"f = {fs[c]:g}", fontsize=10)
            if c == 0:
                ax.set_ylabel(LABELS[r] + "\n dX/dt", fontsize=8)
            if r == 2:
                ax.set_xlabel("X", fontsize=9)
    fig.suptitle(title, fontsize=11)
    p = f"{OUT}/{fname}"; fig.savefig(p, dpi=130); plt.close(fig); return p


def _fli_grid(D, mapk, xsk, vsk, fk, title, fname):
    fs = g(D["monostable"], fk).ravel()
    fig, axes = plt.subplots(3, len(fs), figsize=(11, 8.2), constrained_layout=True)
    im = None
    for r, nm in enumerate(NAMES):
        C = D[nm]; xs = g(C, xsk).ravel(); vs = g(C, vsk).ravel(); eta = gs(C, "eta"); Um = gs(C, "Umin")
        maps = g(C, mapk); XG, VG = np.meshgrid(xs, vs); H0 = .5 * VG ** 2 + U(XG, eta) - Um
        for c in range(len(fs)):
            ax = axes[r, c]
            im = ax.imshow(np.clip(maps[c], 0, 20), origin="lower", aspect="auto",
                           extent=[xs[0], xs[-1], vs[0], vs[-1]], cmap=FLI, vmin=0, vmax=20)
            separatrix(ax, C, None, None); ax.contour(XG, VG, H0, [0.5], colors="w", linewidths=.4, linestyles=":")
            frac = float(np.mean(maps[c][H0 <= 0.5] <= 8))
            ax.text(.03, .96, f"regular {100*frac:.0f}%", transform=ax.transAxes, va="top", color="w", fontsize=8)
            if r == 0:
                ax.set_title(f"f = {fs[c]:g}", fontsize=10)
            if c == 0:
                ax.set_ylabel(LABELS[r] + "\n dX₀/dt", fontsize=8)
            if r == 2:
                ax.set_xlabel("X₀", fontsize=9)
    fig.colorbar(im, ax=axes, shrink=.5, label="FLI após 400 períodos (escuro: toros KAM; claro: caos)")
    fig.suptitle(title, fontsize=12)
    p = f"{OUT}/{fname}"; fig.savefig(p, dpi=130); plt.close(fig); return p


def fig_K2(D):
    return _sections(D, "sec_pts", "sec_fli", "sec_E0", "sec_f",
                     "K2 — Quebra dos toros: seções de Poincaré (1.5 GL, absorvedor congelado)", "K2_secoes_preview.png")


def fig_K3(D):
    return _fli_grid(D, "fli_maps", "fli_xs", "fli_vs", "fli_f",
                     "K3 — Mapas FLI sobre as condições iniciais (1.5 GL)", "K3_FLI_preview.png")


def fig_K4(D):
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    cols = [(.11, .49, .84), (.91, .35, .05), (.18, .62, .27)]
    for k, nm in enumerate(NAMES):
        C = D[nm]; f = g(C, "frac_f").ravel()
        ax.fill_between(f, g(C, "frac_fli6").ravel(), g(C, "frac_fli10").ravel(), color=cols[k], alpha=.18)
        ax.plot(f, g(C, "frac_regular").ravel(), "o-", color=cols[k], ms=4, lw=1.3, label=LABELS[k])
    C = D["shallow_wells"]
    ax.plot(g(C, "abs_f").ravel(), g(C, "abs_regular").ravel(), "s--", color=(.48, .17, .75), ms=6,
            mfc="w", lw=1, label="shallow + absorvedor (2.5 GL)")
    ax.set_xscale("log"); ax.set_xlim(8e-4, .35); ax.set_ylim(0, 1.02)
    ax.set_xlabel("amplitude de forçamento f"); ax.set_ylabel("fração regular (toros KAM)")
    ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("K4 — Fração regular do espaço de fase vs f (banda: limiar FLI de 6 a 10)", fontsize=11)
    p = f"{OUT}/K4_fracao_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


def fig_K5(D):
    C = D["shallow_wells"]; af = g(C, "abs_f").ravel(); ff = g(C, "fli_f").ravel()
    fig, axes = plt.subplots(2, len(af), figsize=(9, 6), constrained_layout=True)
    rows = ["1.5 GL (absorvedor congelado)", "2.5 GL (com absorvedor)"]
    im = None
    for c, f in enumerate(af):
        k = int(np.argmin(np.abs(ff - f)))
        for r in range(2):
            ax = axes[r, c]
            if r == 0:
                M = g(C, "fli_maps")[k]; xs = g(C, "fli_xs").ravel(); vs = g(C, "fli_vs").ravel()
            else:
                M = g(C, "abs_maps")[c]; xs = g(C, "abs_xs").ravel(); vs = g(C, "abs_vs").ravel()
            eta = gs(C, "eta"); Um = gs(C, "Umin"); XG, VG = np.meshgrid(xs, vs); H0 = .5 * VG ** 2 + U(XG, eta) - Um
            im = ax.imshow(np.clip(M, 0, 20), origin="lower", aspect="auto",
                           extent=[xs[0], xs[-1], vs[0], vs[-1]], cmap=FLI, vmin=0, vmax=20)
            separatrix(ax, C, None, None); ax.contour(XG, VG, H0, [0.5], colors="w", linewidths=.4, linestyles=":")
            frac = float(np.mean(M[H0 <= 0.5] <= 8))
            ax.text(.03, .96, f"regular {100*frac:.0f}%", transform=ax.transAxes, va="top", color="w", fontsize=8)
            if r == 0:
                ax.set_title(f"f = {f:g}", fontsize=10)
            if c == 0:
                ax.set_ylabel(rows[r] + "\n dX₀/dt", fontsize=8)
            if r == 1:
                ax.set_xlabel("X₀", fontsize=9)
    fig.colorbar(im, ax=axes, shrink=.5, label="FLI após 400 períodos")
    fig.suptitle("K5 — Acoplamento ao absorvedor encolhe a região KAM (shallow wells)", fontsize=11)
    p = f"{OUT}/K5_absorvedor_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


def fig_K6(D):
    rows = ["shallow_wells", "deep_wells"]; rlab = [LABELS[1], LABELS[2]]
    z1 = g(D["shallow_wells"], "damp_zeta1").ravel(); n = len(z1)
    fig, axes = plt.subplots(2, n, figsize=(10, 5.6), constrained_layout=True)
    for r, nm in enumerate(rows):
        C = D[nm]; xw = gs(C, "xw"); vw = gs(C, "vw")
        early = g(C, "damp_early"); late = g(C, "damp_late"); E0 = g(C, "damp_E0").ravel(); zz = g(C, "damp_zeta1").ravel()
        for c in range(n):
            ax = axes[r, c]; e = early[c]; l = late[c]
            ax.scatter(e[:, 0, :].ravel(), e[:, 1, :].ravel(), s=.2, c=[[.81, .81, .81]], linewidths=0)
            cc = np.repeat(np.clip(E0 / .8, 0, 1)[None, :], l.shape[0], axis=0)
            sz = .3 if zz[c] == 0 else 3
            ax.scatter(l[:, 0, :].ravel(), l[:, 1, :].ravel(), s=sz, c=cc.ravel(), cmap=ENE, vmin=0, vmax=1, linewidths=0)
            separatrix(ax, C, None, None); ax.set_xlim(-xw, xw); ax.set_ylim(-vw, vw)
            if r == 0:
                ax.set_title("ζ₁=0" if zz[c] == 0 else f"ζ₁=1e{int(round(np.log10(zz[c])))}", fontsize=9)
            if c == 0:
                ax.set_ylabel(rlab[r] + "\n dX/dt", fontsize=8)
            if r == 1:
                ax.set_xlabel("X", fontsize=9)
    fig.suptitle("K6 — Amortecimento destrói os toros KAM (1.5 GL): colapso em atratores", fontsize=11)
    p = f"{OUT}/K6_amortecimento_preview.png"; fig.savefig(p, dpi=130); plt.close(fig); return p


def fig_CMP(D):
    """Barras: fração regular 1.5 GL vs 2.5 GL, por caso e forçamento (para o relatório)."""
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), constrained_layout=True, sharey=True)
    for j, nm in enumerate(NAMES):
        C = D[nm]; eta = gs(C, "eta"); Um = gs(C, "Umin")
        ff = g(C, "fli_f").ravel(); af = g(C, "abs_f").ravel()
        fmaps = g(C, "fli_maps"); fxs = g(C, "fli_xs").ravel(); fvs = g(C, "fli_vs").ravel()
        XG, VG = np.meshgrid(fxs, fvs); H0 = .5 * VG ** 2 + U(XG, eta) - Um
        areg = g(C, "abs_regular").ravel()
        r15 = [100 * float(np.mean(fmaps[int(np.argmin(np.abs(ff - f)))][H0 <= .5] <= 8)) for f in af]
        r25 = [100 * v for v in areg]
        x = np.arange(len(af)); w = .38; ax = axes[j]
        b1 = ax.bar(x - w / 2, r15, w, label="1.5 GL (QZS)", color="#1f4e8c")
        b2 = ax.bar(x + w / 2, r25, w, label="2.5 GL (QZS-ADV)", color="#bc4b51")
        for b in list(b1) + list(b2):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1, f"{b.get_height():.0f}",
                    ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels([f"f={f:g}" for f in af]); ax.set_ylim(0, 108)
        ax.set_title(LABELS[j].split("  ")[0]); ax.grid(axis="y", alpha=.3)
        if j == 0:
            ax.set_ylabel("fração regular (%) — toros KAM"); ax.legend(loc="lower left", fontsize=8)
    fig.suptitle("Fração regular do espaço de fase: QZS (1.5 GL) vs QZS-ADV (2.5 GL)", fontsize=12)
    p = f"{OUT}/CMP_fracao_regular_preview.png"; fig.savefig(p, dpi=150); plt.close(fig); return p


if __name__ == "__main__":
    D = load()
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else None
    allfigs = {"K1": fig_K1, "K2": fig_K2, "K3": fig_K3, "K4": fig_K4, "K5": fig_K5,
               "K6": fig_K6, "K7": fig_K7, "K8": fig_K8, "K9": fig_K9, "CMP": fig_CMP}
    for key, fn in allfigs.items():
        if only and key not in only:
            continue
        print(fn(D), flush=True)
