"""Como a sintonia β ótima muda com a profundidade do poço.

Sobrepõe as curvas fração-regular-KAM(β) dos três casos (de beta_tuning.py, Ω=1) e
sombreia a banda de frequência ω(E) do núcleo (dE≤ECUT) de cada caso — a "zona de
perigo" onde ω(E)=β destrói os toros. Mostra que poços rasos têm um vale de caos largo
em β baixo-médio, enquanto poços profundos são robustos exceto por um pico estreito em
β alto (onde ω(E)=β encontra o núcleo profundo).

Requer beta_tuning_result.json (rode beta_tuning.py --full antes) e dados_kam.mat.
"""
import json

import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAT = "/home/user/Chaos/dados_kam.mat"
RES = "/home/user/Chaos/work/beta_tuning_result.json"
OUT = "/home/user/Chaos/work/beta_vs_depth_preview.png"
ECUT = 0.5
COLS = {"monostable": "#1f4e8c", "shallow_wells": "#e0a020", "deep_wells": "#bc4b51"}
LAB = {"monostable": "Monostable (ω₀=0)", "shallow_wells": "Shallow wells (ω₀=0,65)", "deep_wells": "Deep wells (ω₀=1,56)"}


def unwrap(x):
    while isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
        x = x.reshape(-1)[0]
    return x


def gs(C, k):
    return float(np.asarray(getattr(C, k)).ravel()[0])


def main():
    D = unwrap(loadmat(MAT, squeeze_me=False, struct_as_record=False)["D"])
    bands = {}
    for nm in COLS:
        C = unwrap(getattr(D, nm)); dE = np.asarray(C.dE).ravel(); w = np.asarray(C.omega_quad).ravel()
        reg = dE <= ECUT
        bands[nm] = (float(w[reg].min()), float(w[reg].max()), gs(C, "omega0_sqrtK"))
    R = json.load(open(RES))
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    for nm in COLS:
        d = R["cases"][nm]; bs = np.array(d["betas"]); reg = 100 * np.array(d["regular"])
        ax.plot(bs, reg, "o-", color=COLS[nm], lw=1.8, ms=3.5, label=LAB[nm])
        lo, hi, _ = bands[nm]
        ax.axvspan(lo, hi, color=COLS[nm], alpha=.06)
    ax.axvline(1.0, color="#1f7a1f", ls="--", lw=1.2)
    ax.text(1.0, 103, "β=Ω=1 (antirress.)", color="#1f7a1f", fontsize=8, ha="center")
    ax.set_xlabel("β (sintonia do absorvedor)"); ax.set_ylabel("fração regular KAM (%)")
    ax.set_ylim(0, 108); ax.set_xlim(0.1, 1.9); ax.grid(alpha=.25); ax.legend(loc="upper center", fontsize=9)
    ax.set_title("Sintonia β vs profundidade do poço (Ω=1, f=0,05); faixas = banda de ω(E) do núcleo")
    fig.savefig(OUT, dpi=150); print("salvo", OUT)
    for nm in COLS:
        lo, hi, w0 = bands[nm]
        print(f"  {nm:14}: ω0={w0:.2f}, ω(E) núcleo ∈ [{lo:.2f},{hi:.2f}]")


if __name__ == "__main__":
    main()
