"""Classificação regular/caótico independente da dimensão, via LEI DE CRESCIMENTO do FLI.

O problema: o FLI de uma órbita regular cresce como log(t); o de uma caótica, como
lambda*t. O VALOR absoluto do FLI regular é maior em 4D (2.5 GL) do que em 2D (1.5 GL)
porque o toro tem mais frequências — então um limiar fixo (8, 10) conta caos a mais no
2.5 GL. A cura: classificar pela taxa de crescimento entre t e 2t, que é ~log10(2)≈0.30
para QUALQUER órbita regular (independente da dimensão) e ordens de grandeza maior para
caótica. Assim a comparação 1.5 GL x 2.5 GL fica no mesmo pé.

Grava o FLI em snapshots (períodos SNAPS) e define:
    caótico  <=>  FLI(2T) - FLI(T) > TAU     (crescimento super-logarítmico)
com T = SNAPS[0], 2T = SNAPS[1]. TAU padrão 1.0 (regular ~0.3; margem folgada).

Uso:
  python recalibrate_fli.py            # preview (grade/períodos pequenos)
  python recalibrate_fli.py --full     # produção (nx=160, 400 períodos)
"""
import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import kam_compute as K

NAMES = ["monostable", "shallow_wells", "deep_wells"]
FS = [0.01, 0.05, 0.15]
TAU = float(os.environ.get("KAM_TAU", "1.0"))
OUT = os.environ.get("KAM_OUT", "/home/user/Chaos/work/recalibrate_result.json")

PRESET = {"preview": dict(nx=48, snaps=(60, 120), ns=40),
          "full":    dict(nx=160, snaps=(200, 400), ns=100)}


def gridpts(name, nx):
    m = K.CASES[name]
    xs = np.linspace(-m["xw"], m["xw"], nx); vs = np.linspace(-m["vw"], m["vw"], nx)
    X, V = np.meshgrid(xs, vs)
    return xs, vs, X.ravel(), V.ravel()


def snap_1dof(x, v, eta, f, snaps, ns):
    """FLI (2D tangente) do modelo 1.5 GL, gravado nos períodos `snaps`."""
    n = x.size; h = K.T / ns; t = 0.
    dx = np.ones(n) / np.sqrt(2); dv = dx.copy(); logg = np.zeros(n)
    out = {}; smax = max(snaps)
    for c in range(smax):
        for s in range(ns):
            for k in range(3):
                x = x + K.CY[k] * h * v; t += K.CY[k] * h; dx = dx + K.CY[k] * h * dv
                v = v + K.DY[k] * h * (f * np.cos(t) - K.dU(x, eta)); dv = dv - K.DY[k] * h * K.d2U(x, eta) * dx
            x = x + K.CY[3] * h * v; t += K.CY[3] * h; dx = dx + K.CY[3] * h * dv
        t = (c + 1) * K.T
        nrm = np.sqrt(dx * dx + dv * dv); logg += np.log10(nrm); dx /= nrm; dv /= nrm
        if (c + 1) in snaps:
            out[c + 1] = logg.copy()
    return out


def snap_2dof(x, v, x2, v2, eta, f, snaps, ns):
    """FLI (4D tangente) do modelo 2.5 GL, gravado nos períodos `snaps`."""
    n = x.size; h = K.T / ns; t = 0.; mb = K.MU * K.BETA2
    d = [np.ones(n) / 2 for _ in range(4)]; logg = np.zeros(n)
    out = {}; smax = max(snaps)
    for c in range(smax):
        for s in range(ns):
            for k in range(4):
                x = x + K.CY[k] * h * v; x2 = x2 + K.CY[k] * h * v2; t += K.CY[k] * h
                d[0] = d[0] + K.CY[k] * h * d[1]; d[2] = d[2] + K.CY[k] * h * d[3]
                if k == 3:
                    break
                coup = mb * (x - x2)
                v = v + K.DY[k] * h * (f * np.cos(t) - K.dU(x, eta) - coup); v2 = v2 + K.DY[k] * h * (K.BETA2 * (x - x2))
                k2 = K.d2U(x, eta)
                d[1] = d[1] + K.DY[k] * h * (-(k2 + mb) * d[0] + mb * d[2]); d[3] = d[3] + K.DY[k] * h * (K.BETA2 * (d[0] - d[2]))
        t = (c + 1) * K.T
        nrm = np.sqrt(d[0]**2 + d[1]**2 + d[2]**2 + d[3]**2); logg += np.log10(nrm); d = [q / nrm for q in d]
        if (c + 1) in snaps:
            out[c + 1] = logg.copy()
    return out


def job(spec):
    model, nm, f, cfg = spec
    eta = K.CASES[nm]["eta"]; xs, vs, X, V = gridpts(nm, cfg["nx"]); snaps = cfg["snaps"]
    if model == "1.5":
        snp = snap_1dof(X.copy(), V.copy(), eta, f, snaps, cfg["ns"])
    else:
        snp = snap_2dof(X.copy(), V.copy(), X.copy(), V.copy(), eta, f, snaps, cfg["ns"])
    H0 = .5 * V**2 + K.U(X, eta) - K.U(K.xmin(eta), eta)
    inside = H0 <= K.ECUT
    t1, t2 = snaps
    growth = snp[t2] - snp[t1]                     # crescimento log-super-linear
    fli_abs = snp[t2]                               # valor absoluto no tempo final
    reg_growth = float(np.mean(growth[inside] <= TAU))       # novo classificador
    reg_abs8 = float(np.mean(fli_abs[inside] <= 8))          # limiar absoluto antigo (8)
    return (model, nm, f, reg_growth, reg_abs8,
            float(np.median(growth[inside])), float(np.percentile(growth[inside], 90)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("KAM_WORKERS", "4")))
    args = ap.parse_args()
    cfg = PRESET["full" if args.full else "preview"]
    mode = "full" if args.full else "preview"
    print(f"[{mode}] TAU={TAU}  snaps={cfg['snaps']}  nx={cfg['nx']}  workers={args.workers}", flush=True)
    specs = [(model, nm, f, cfg) for model in ("1.5", "2.5") for nm in NAMES for f in FS]
    res = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [pool.submit(job, s) for s in specs]
        for k, fu in enumerate(as_completed(futs), 1):
            r = fu.result(); res[(r[0], r[1], r[2])] = r
            print(f"  {k}/{len(specs)} {r[0]}GL {r[1]} f={r[2]}: reg(cresc)={100*r[3]:.0f}% reg(abs8)={100*r[4]:.0f}%", flush=True)
    print("\n===== COMPARAÇÃO: fração regular, limiar ABSOLUTO(8) vs LEI DE CRESCIMENTO =====")
    print(f"{'caso':14}{'f':>6}  | {'1.5GL abs':>9} {'2.5GL abs':>9} {'Δabs':>6}  | {'1.5GL cre':>9} {'2.5GL cre':>9} {'Δcre':>6}")
    table = []
    for nm in NAMES:
        for f in FS:
            a1 = res[("1.5", nm, f)][4]; a2 = res[("2.5", nm, f)][4]
            g1 = res[("1.5", nm, f)][3]; g2 = res[("2.5", nm, f)][3]
            print(f"{nm:14}{f:>6g}  | {100*a1:>8.0f}% {100*a2:>8.0f}% {100*(a2-a1):>+5.0f}  | {100*g1:>8.0f}% {100*g2:>8.0f}% {100*(g2-g1):>+5.0f}")
            table.append(dict(case=nm, f=f, abs_1p5=a1, abs_2p5=a2, grow_1p5=g1, grow_2p5=g2))
    json.dump(dict(tau=TAU, snaps=list(cfg["snaps"]), mode=mode, rows=table), open(OUT, "w"), indent=1)
    print(f"\nsalvo {OUT}")


if __name__ == "__main__":
    main()
