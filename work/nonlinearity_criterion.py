"""Critério de não linearidade x diferença de caos entre QZS (1.5 GL) e QZS-ADV (2.5 GL).

Usa a classificação de caos por LEI DE CRESCIMENTO do FLI (recalibrate_fli.py) como
métrica de caos, e liga a DIFERENÇA de caos entre os dois sistemas à não linearidade
do primário, quantificada pela relação frequência-amplitude ω(E).

Tese: o caos extra do 2.5 GL vem da ressonância interna primário↔absorvedor, quando
a frequência do primário ω(E) iguala a frequência natural do absorvedor ω_a = β. O
QZS achata o potencial → rigidez local ω0 baixa → a frequência de operação (baixa
energia, onde vivem as órbitas regulares) já está próxima de ω_a → a ressonância cai
DENTRO do núcleo regular e o injeta caos. Quanto MENOR ω0 (mais forte a não linearidade
do QZS no ponto de operação), MAIOR o aumento de caos ao acoplar o absorvedor.

Critério operacional: a energia dE_res do cruzamento ω(E) = ω_a prediz a amplificação
do caos — dE_res no núcleo (baixa energia) ⇒ grande amplificação; dE_res na barreira
⇒ pequena.

Gera nonlinearity_criterion.png e imprime a tabela do critério.
"""
import json

import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAT = "/home/user/Chaos/dados_kam.mat"
RECAL = "/home/user/Chaos/work/recalibrate_result.json"
OUT = "/home/user/Chaos/work/nonlinearity_criterion_preview.png"
BETA = float(np.sqrt(0.1225))   # ω_a: frequência natural do absorvedor
ECUT = 0.5
NAMES = ["monostable", "shallow_wells", "deep_wells"]
LAB = {"monostable": "Monostable", "shallow_wells": "Shallow wells", "deep_wells": "Deep wells"}
COL = {"monostable": "#1f4e8c", "shallow_wells": "#e0a020", "deep_wells": "#bc4b51"}


def unwrap(x):
    while isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
        x = x.reshape(-1)[0]
    return x


def g(C, k):
    return np.asarray(getattr(C, k))


def gs(C, k):
    return float(g(C, k).ravel()[0])


def resonance_dE(dE, w, over, target):
    """Menor energia dE onde ω(E)=target (cruzamento interpolado)."""
    best = None
    for b in np.unique(over):
        m = over == b
        o = np.argsort(dE[m]); de = dE[m][o]; ww = w[m][o]
        s = ww - target
        for i in np.where(s[:-1] * s[1:] < 0)[0]:
            dc = de[i] + (de[i + 1] - de[i]) * (-s[i]) / (s[i + 1] - s[i])
            if best is None or dc < best:
                best = float(dc)
    return best


def main():
    D = unwrap(loadmat(MAT, squeeze_me=False, struct_as_record=False)["D"])
    recal = {(r["case"], r["f"]): r for r in json.load(open(RECAL))["rows"]}

    info = {}
    for nm in NAMES:
        C = unwrap(getattr(D, nm))
        dE = g(C, "dE").ravel(); w = g(C, "omega_quad").ravel(); over = g(C, "over").ravel()
        w0 = gs(C, "omega0_sqrtK")
        dres = resonance_dE(dE, w, over, BETA)
        # caos = 1 - fração regular (lei de crescimento); aumento ao adicionar o absorvedor
        dchaos = {f: 100 * (recal[(nm, f)]["grow_1p5"] - recal[(nm, f)]["grow_2p5"]) for f in (0.01, 0.05, 0.15)}
        info[nm] = dict(dE=dE, w=w, over=over, w0=w0, dres=dres, dchaos=dchaos)

    # ---- figura ----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # painel esquerdo: ω(E) com ω_a e a ressonância interna
    axL.axhspan(BETA * 0.9, BETA * 1.1, color="0.85", zorder=0)
    axL.axhline(BETA, color="0.4", ls="--", lw=1)
    axL.text(1.3, BETA, "  ω_a = β (absorvedor)", va="center", fontsize=9, color="0.3")
    axL.axvline(ECUT, color="0.6", ls=":", lw=.8)
    axL.text(ECUT, 1.72, " núcleo regular (dE ≤ 0,5) →", fontsize=8, color="0.4", ha="right", rotation=90, va="top")
    for nm in NAMES:
        d = info[nm]
        for b in np.unique(d["over"]):
            m = d["over"] == b; o = np.argsort(d["dE"][m])
            axL.plot(d["dE"][m][o], d["w"][m][o], color=COL[nm], lw=1.6,
                     label=LAB[nm] if b == np.unique(d["over"])[0] else None)
        if d["dres"] is not None:
            axL.plot([d["dres"]], [BETA], "o", ms=9, mfc=COL[nm], mec="w", mew=1.2, zorder=5)
    axL.set_xscale("log"); axL.set_xlim(1e-4, 1.7); axL.set_ylim(0, 1.75)
    axL.set_xlabel("E − U_min  (energia acima do mínimo)"); axL.set_ylabel("ω(E) / Ω  (frequência do primário)")
    axL.set_title("A) Onde a frequência do primário cruza a do absorvedor")
    axL.legend(loc="upper left", fontsize=9)

    # painel direito: critério — dE_res vs aumento de caos (lei de crescimento)
    for nm in NAMES:
        d = info[nm]
        axR.scatter([d["dres"]], [d["dchaos"][0.05]], s=160, color=COL[nm], edgecolor="w", zorder=5)
        axR.annotate(f"{LAB[nm]}\n(ω₀={d['w0']:.2f})", (d["dres"], d["dchaos"][0.05]),
                     textcoords="offset points", xytext=(12, -4), fontsize=9)
    axR.set_xscale("log"); axR.set_xlabel("energia da ressonância interna  dE_res  (ω(E)=ω_a)")
    axR.set_ylabel("aumento de caos 1.5→2.5 GL  (p.p., f=0,05)")
    axR.set_title("B) Critério: ressonância no núcleo ⇒ muito caos extra")
    axR.grid(alpha=.3); axR.set_ylim(-5, 60)
    axR.axvspan(1e-4, 0.05, color="#bc4b51", alpha=.06)
    axR.axvspan(0.05, 1.0, color="#1f7a1f", alpha=.06)
    axR.text(0.006, 5, "ressonância\nno núcleo", fontsize=8, color="#bc4b51", ha="center")
    axR.text(0.3, 45, "ressonância\nna barreira", fontsize=8, color="#1f7a1f", ha="center")

    fig.suptitle("Critério não linearidade × diferença de caos: QZS (1.5 GL) vs QZS-ADV (2.5 GL)", fontsize=12)
    fig.savefig(OUT, dpi=150)
    print("salvo", OUT)

    print("\n===== CRITÉRIO =====")
    print(f"ω_a (absorvedor) = {BETA:.3f}\n")
    print(f"{'caso':14}{'ω0':>7}{'dE_res':>9}{'  aumento de caos (p.p.)':>26}")
    for nm in NAMES:
        d = info[nm]
        dc = "  ".join(f"f={f}:{d['dchaos'][f]:+.0f}" for f in (0.01, 0.05, 0.15))
        print(f"{nm:14}{d['w0']:>7.2f}{d['dres']:>9.4f}   {dc}")


if __name__ == "__main__":
    main()
