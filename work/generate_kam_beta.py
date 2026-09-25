"""Mapas KAM (FLI 2.5 GL conservativo) para diferentes sintonias β do absorvedor.

Mostra visualmente como os toros KAM mudam com a sintonia: β=0.35 (baseline),
β=0.70 (vale de caos) e β=1.00 (=Ω, antirressonância). Linhas = casos, colunas = β.
Cor = FLI absoluto após 400 períodos (como nas figuras K7); anota a fração regular
pela lei de crescimento (dimensão-independente).
"""
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import kam_compute as K

NAMES = ["monostable", "shallow_wells", "deep_wells"]
LAB = {"monostable": "Monostable  (ω₀=0)", "shallow_wells": "Shallow wells  (ω₀=0,65)", "deep_wells": "Deep wells  (ω₀=1,56)"}
BETAS = [(0.35, "β=0,35 (baseline)"), (0.70, "β=0,70 (vale de caos)"), (1.00, "β=1,00  (β=Ω, antirressonância)")]
FFORCE = 0.05
NX = 120
SNAPS = (200, 400)
NS = 100
FLI = LinearSegmentedColormap.from_list("f", ["#" + h for h in
      ["0b1d3a", "1f4e8c", "3a8fb7", "9ad0c2", "f4e285", "f4a259", "bc4b51", "5b1a18"]])


def snap_map(nm, beta):
    """FLI 4D conservativo (β dado); retorna xs, vs, mapa FLI(400 per) e fração regular (lei de cresc.)."""
    eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; b2 = beta ** 2; mb = K.MU * b2
    xs = np.linspace(-m["xw"], m["xw"], NX); vs = np.linspace(-m["vw"], m["vw"], NX)
    X, V = np.meshgrid(xs, vs); x = X.ravel(); v = V.ravel(); x2 = x.copy(); v2 = v.copy()
    n = x.size; h = K.T / NS; t = 0.
    d = [np.ones(n) / 2 for _ in range(4)]; logg = np.zeros(n); snaps = {}
    for c in range(max(SNAPS)):
        for s in range(NS):
            for k in range(4):
                x = x + K.CY[k] * h * v; x2 = x2 + K.CY[k] * h * v2; t += K.CY[k] * h
                d[0] = d[0] + K.CY[k] * h * d[1]; d[2] = d[2] + K.CY[k] * h * d[3]
                if k == 3:
                    break
                coup = mb * (x - x2)
                v = v + K.DY[k] * h * (FFORCE * np.cos(t) - K.dU(x, eta) - coup); v2 = v2 + K.DY[k] * h * (b2 * (x - x2))
                k2 = K.d2U(x, eta)
                d[1] = d[1] + K.DY[k] * h * (-(k2 + mb) * d[0] + mb * d[2]); d[3] = d[3] + K.DY[k] * h * (b2 * (d[0] - d[2]))
        t = (c + 1) * K.T
        nrm = np.sqrt(d[0]**2 + d[1]**2 + d[2]**2 + d[3]**2); logg += np.log10(nrm); d = [q / nrm for q in d]
        if (c + 1) in SNAPS:
            snaps[c + 1] = logg.copy()
    H0 = .5 * V.ravel()**2 + K.U(x, eta) - K.U(K.xmin(eta), eta)
    inside = H0 <= K.ECUT
    growth = snaps[SNAPS[1]] - snaps[SNAPS[0]]
    reg = float(np.mean(growth[inside] <= 1.0))
    return nm, beta, xs, vs, snaps[SNAPS[1]].reshape(NX, NX), reg


def main():
    jobs = [(nm, b) for nm in NAMES for b, _ in BETAS]
    out = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("KAM_WORKERS", "4"))) as pool:
        futs = [pool.submit(snap_map, nm, b) for nm, b in jobs]
        for k, fu in enumerate(as_completed(futs), 1):
            nm, beta, xs, vs, mp, reg = fu.result()
            out[(nm, beta)] = (xs, vs, mp, reg)
            print(f"  {k}/{len(jobs)} {nm} β={beta}: regular(cresc)={100*reg:.0f}%", flush=True)

    fig, axes = plt.subplots(3, 3, figsize=(11, 10.5), constrained_layout=True)
    im = None
    for r, nm in enumerate(NAMES):
        eta = K.CASES[nm]["eta"]; Um = K.U(K.xmin(eta), eta); depth = K.U(0., eta) - Um
        for c, (beta, blab) in enumerate(BETAS):
            xs, vs, mp, reg = out[(nm, beta)]
            ax = axes[r, c]
            im = ax.imshow(np.clip(mp, 0, 20), origin="lower", aspect="auto",
                           extent=[xs[0], xs[-1], vs[0], vs[-1]], cmap=FLI, vmin=0, vmax=20)
            XG, VG = np.meshgrid(xs, vs); HG = .5 * VG**2 + K.U(XG, eta) - Um
            if depth > 0:
                ax.contour(XG, VG, HG, [depth], colors="k", linewidths=.5, linestyles="--")
            ax.contour(XG, VG, HG, [0.5], colors="w", linewidths=.4, linestyles=":")
            ax.text(.03, .96, f"regular {100*reg:.0f}%", transform=ax.transAxes, va="top", color="w", fontsize=9)
            if r == 0:
                ax.set_title(blab, fontsize=10)
            if c == 0:
                ax.set_ylabel(LAB[nm] + "\n dX₀/dt", fontsize=9)
            if r == 2:
                ax.set_xlabel("X₀", fontsize=9)
    fig.colorbar(im, ax=axes, shrink=.4, label="FLI após 400 períodos (escuro: toros KAM; claro: caos)")
    fig.suptitle("Toros KAM (2.5 GL conservativo, f=0,05) vs sintonia β do absorvedor", fontsize=13)
    p = "/home/user/Chaos/work/KAM_beta_tuning_preview.png"; fig.savefig(p, dpi=140); print("salvo", p)


if __name__ == "__main__":
    main()
