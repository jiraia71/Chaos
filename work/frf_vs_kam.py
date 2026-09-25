"""Antirressonância (FRF linearizado) x preservação dos toros KAM, em função de β.

Junta duas curvas vs β do absorvedor:
  - amplitude do primário no forçamento (FRF linearizado, Ω=1): mínimo = antirressonância;
  - fração regular KAM do sistema conservativo 2.5 GL (de beta_tuning.py).
Mostra que a antirressonância (β≈Ω=1) coincide com a janela em que o KAM se recupera,
enquanto o baseline β=0.35 e o vale β≈0.5–0.9 são ruins nos dois critérios.

Requer beta_tuning_result.json (rode beta_tuning.py --full antes).
"""
import json

import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MU, W, Z1, Z2 = 0.1, 1.0, 0.01, 0.01   # μ, Ω, ζ1, ζ2
MAT = "/home/user/Chaos/dados_kam.mat"
RES = "/home/user/Chaos/work/beta_tuning_result.json"
OUT = "/home/user/Chaos/work/FRF_vs_KAM_preview.png"
NAMES = ["monostable", "shallow_wells", "deep_wells"]
LAB = {"monostable": "Monostable", "shallow_wells": "Shallow wells", "deep_wells": "Deep wells"}


def unwrap(x):
    while isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
        x = x.reshape(-1)[0]
    return x


def primary_amp(w0, beta):
    """|X/f| do primário linearizado com absorvedor, no forçamento Ω."""
    b2 = beta ** 2
    den = (w0**2 - W**2 + 2j*Z1*W + MU*b2) - MU*b2*b2 / (b2 - W**2 + 2j*Z2*W)
    return 1.0 / abs(den)


def main():
    D = unwrap(loadmat(MAT, squeeze_me=False, struct_as_record=False)["D"])
    R = json.load(open(RES))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), constrained_layout=True)
    for j, nm in enumerate(NAMES):
        C = unwrap(getattr(D, nm)); w0 = float(np.asarray(getattr(C, "omega0_sqrtK")).ravel()[0])
        d = R["cases"][nm]; bs = np.array(d["betas"]); reg = 100 * np.array(d["regular"])
        ax = axes[j]
        ax.plot(bs, reg, "o-", color="#1f4e8c", lw=1.8, ms=4)
        ax.set_ylim(0, 100); ax.set_xlabel("β (sintonia do absorvedor)")
        if j == 0:
            ax.set_ylabel("fração regular KAM (%)", color="#1f4e8c")
        ax.tick_params(axis="y", labelcolor="#1f4e8c")
        axr = ax.twinx(); bb = np.linspace(bs.min(), bs.max(), 400)
        amp = np.array([primary_amp(w0, b) for b in bb])
        axr.plot(bb, amp, color="#bc4b51", lw=1.8); axr.set_ylim(0, max(amp) * 1.05)
        if j == 2:
            axr.set_ylabel("amplitude primário @Ω=1  (↓ = mais antirressonância)", color="#bc4b51")
        axr.tick_params(axis="y", labelcolor="#bc4b51")
        ax.axvline(0.35, color="0.5", ls=":", lw=1); ax.text(0.35, 95, "β=0.35\n(baseline)", fontsize=7, ha="center", color="0.4")
        ax.axvline(1.0, color="#1f7a1f", ls="--", lw=1.2); ax.text(1.0, 95, "β≈1\n(antirress.)", fontsize=7, ha="center", color="#1f7a1f")
        low = reg < 10
        if low.any():
            ax.axvspan(bs[low].min(), bs[low].max(), color="#bc4b51", alpha=.08)
        ax.set_title(f"{LAB[nm]}  (ω₀={w0:.2f})"); ax.grid(alpha=.25)
    fig.suptitle("Sintonia β: antirressonância (FRF) × preservação dos toros KAM  —  f=0.05, ζ=0.01, Ω=1", fontsize=12)
    fig.savefig(OUT, dpi=150); print("salvo", OUT)


if __name__ == "__main__":
    main()
