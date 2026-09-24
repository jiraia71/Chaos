"""Augment an existing dados_kam.mat with the full 2.5-DOF (QZS-ADV) analysis.

The canonical pipeline is kam_compute.py (compute the cache) + export_kam_matlab.py
(assemble the .mat). That recomputes everything, including the already-validated
1.5-DOF stages. This script is the fast path: it keeps every field of an existing
dados_kam.mat verbatim and only *adds* the new 2.5-DOF fields, computed here in
memory by reusing kam_compute's integrators:

  abs_*    FLI maps with the absorber, for all three cases (K7)
  sec2_*   stroboscopic Poincare sections of the 2.5-DOF system (K8)
  damp2_*  fate of the tori with damping (zeta1, zeta2) in 2.5-DOF (K9)

Usage:
  python augment_adv.py            # preview: coarse grids/short runs (fast smoke test)
  python augment_adv.py --full     # production resolution (heavy; use many workers)

Env: KAM_SRC / KAM_OUT override the input/output .mat paths.
"""
import argparse
import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

import kam_compute as K

SRC = Path(os.environ.get("KAM_SRC", "/home/user/Chaos/dados_kam.mat"))

# resolution knobs: (map grid nx, map periods, map steps/period,
#                    section periods, section steps, damping periods, damping steps)
PRESET = {
    "preview": dict(nx=48, abs_p=60, abs_ns=40, sec_p=200, sec_ns=120, dmp_p=300, dmp_ns=180),
    "full":    dict(nx=160, abs_p=400, abs_ns=100, sec_p=1500, sec_ns=200, dmp_p=3000, dmp_ns=180),
}


def unwrap(x):
    """Nested MATLAB structs come back wrapped in object ndarrays; peel to the
    mat_struct so getattr(...) works."""
    while isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
        x = x.reshape(-1)[0]
    return x


def grid(name, nx):
    m = K.CASES[name]
    xs = np.linspace(-m["xw"], m["xw"], nx)
    vs = np.linspace(-m["vw"], m["vw"], nx)
    X, V = np.meshgrid(xs, vs)
    return xs, vs, X.ravel(), V.ravel()


def compute_abs(name, f, cfg):
    """2.5-DOF FLI map (absorber at rest relative to the primary, t = 0)."""
    eta = K.CASES[name]["eta"]
    xs, vs, X, V = grid(name, cfg["nx"])
    fli, _ = K.yoshida_2dof(X.copy(), V.copy(), X.copy(), V.copy(), eta, f,
                            cfg["abs_p"], cfg["abs_ns"], tag=f"abs {name} f={f}")
    H0 = .5 * V**2 + K.U(X, eta) - K.U(K.xmin(eta), eta)
    nx = cfg["nx"]
    reg = float(np.mean(fli[H0 <= K.ECUT] <= 8))
    return xs, vs, fli.reshape(nx, nx).astype(np.float32), reg


def compute_sec2(name, f, cfg):
    """2.5-DOF stroboscopic section, projected on the primary (X, dX/dt)."""
    eta = K.CASES[name]["eta"]
    x0, v0 = K.ics_line(name)
    fli, pts = K.yoshida_2dof(x0.copy(), v0.copy(), x0.copy(), v0.copy(), eta, f,
                              cfg["sec_p"], cfg["sec_ns"], record=True, tag=f"sec2 {name} f={f}")
    E0 = .5 * v0**2 + K.U(x0, eta) - K.U(K.xmin(eta), eta)
    return pts.astype(np.float32), fli, E0


def compute_damp2(name, f, z1, z2, cfg):
    """2.5-DOF RK4 with viscous damping (zeta1 primary, zeta2 absorber)."""
    eta = K.CASES[name]["eta"]
    x, v = K.ics_line(name)
    x2 = x.copy(); v2 = v.copy()
    ns = cfg["dmp_ns"]; h = K.T / ns; mb = K.MU * K.BETA2
    cs = np.cos(np.arange(2 * ns + 1) * np.pi / ns)
    periods = cfg["dmp_p"]
    pts = np.empty((periods, 2, x.size), np.float32)

    def deriv(x, v, x2, v2, force):
        return (v, force - K.dU(x, eta) - mb * (x - x2) - 2 * z1 * v,
                v2, K.BETA2 * (x - x2) - 2 * z2 * v2)

    for c in range(periods):
        for s in range(ns):
            k1 = deriv(x, v, x2, v2, f * cs[2 * s])
            k2 = deriv(x + .5 * h * k1[0], v + .5 * h * k1[1], x2 + .5 * h * k1[2], v2 + .5 * h * k1[3], f * cs[2 * s + 1])
            k3 = deriv(x + .5 * h * k2[0], v + .5 * h * k2[1], x2 + .5 * h * k2[2], v2 + .5 * h * k2[3], f * cs[2 * s + 1])
            k4 = deriv(x + h * k3[0], v + h * k3[1], x2 + h * k3[2], v2 + h * k3[3], f * cs[2 * s + 2])
            x = x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v = v + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            x2 = x2 + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
            v2 = v2 + h / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        pts[c, 0] = x; pts[c, 1] = v
    E0 = K.U(K.ics_line(name)[0], eta) - K.U(K.xmin(eta), eta)
    return pts, E0


def job(spec):
    """Run one unit of work in a worker process. spec = (kind, nm, idx, cfg[, args])."""
    kind, nm, idx, cfg = spec[0], spec[1], spec[2], spec[3]
    if kind == "abs":
        f = spec[4]
        xs, vs, mp, reg = compute_abs(nm, f, cfg)
        return kind, nm, idx, (xs, vs, mp, reg)
    if kind == "sec2":
        f = spec[4]
        pts, fli, E0 = compute_sec2(nm, f, cfg)
        return kind, nm, idx, (pts, fli, E0)
    if kind == "damp2":
        z1, z2 = spec[4]
        pts, E0 = compute_damp2(nm, .05, z1, z2, cfg)
        return kind, nm, idx, (pts[:150], pts[-300:], E0)
    raise ValueError(kind)


def main():
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="production resolution (heavy)")
    ap.add_argument("--out", default=os.environ.get("KAM_OUT", ""))
    ap.add_argument("--workers", type=int, default=int(os.environ.get("KAM_WORKERS", "4")))
    args = ap.parse_args()
    cfg = PRESET["full" if args.full else "preview"]
    mode = "full" if args.full else "preview"
    out = Path(args.out) if args.out else SRC.with_name(
        "dados_kam.mat" if args.full else "dados_kam_preview.mat")

    print(f"[{mode}] lendo {SRC} (workers={args.workers})", flush=True)
    raw = loadmat(SRC, squeeze_me=False, struct_as_record=False)
    Dsrc = unwrap(raw["D"])
    cases = list(Dsrc._fieldnames)

    # copia verbatim todos os campos 1.5 GL existentes (formas preservadas)
    D = {}
    for nm in cases:
        Csrc = unwrap(getattr(Dsrc, nm))
        D[nm] = {f: getattr(Csrc, f) for f in Csrc._fieldnames}

    # monta a lista de tarefas 2.5 GL e distribui num pool
    specs = []
    for nm in cases:
        for i, f in enumerate(K.F_ABS_ALL):
            specs.append(("abs", nm, i, cfg, f))
        for i, f in enumerate(K.F_SECTIONS):
            specs.append(("sec2", nm, i, cfg, f))
        for i, zz in enumerate(K.DAMP2_ZETA):
            specs.append(("damp2", nm, i, cfg, zz))
    # buffers indexados por caso
    buf = {nm: dict(abs=[None] * len(K.F_ABS_ALL), sec2=[None] * len(K.F_SECTIONS),
                    damp2=[None] * len(K.DAMP2_ZETA)) for nm in cases}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [pool.submit(job, s) for s in specs]
        for k, fu in enumerate(as_completed(futs), 1):
            kind, nm, idx, payload = fu.result()
            buf[nm][kind][idx] = payload
            print(f"  done {k}/{len(specs)}  {kind} {nm}[{idx}]  [{time.time()-t0:.0f}s]", flush=True)

    for nm in cases:
        c = D[nm]
        # abs_* (K7)
        xs = buf[nm]["abs"][0][0]; vs = buf[nm]["abs"][0][1]
        c["abs_f"] = np.array(K.F_ABS_ALL, float).reshape(1, -1)
        c["abs_maps"] = np.stack([p[2] for p in buf[nm]["abs"]]).astype(np.float32)
        c["abs_xs"] = xs.reshape(1, -1); c["abs_vs"] = vs.reshape(-1, 1)
        c["abs_regular"] = np.array([p[3] for p in buf[nm]["abs"]], float).reshape(1, -1)
        # sec2_* (K8)
        c["sec2_f"] = np.array(K.F_SECTIONS, float).reshape(1, -1)
        c["sec2_pts"] = np.stack([p[0] for p in buf[nm]["sec2"]]).astype(np.float32)
        c["sec2_fli"] = np.stack([p[1] for p in buf[nm]["sec2"]])
        c["sec2_E0"] = buf[nm]["sec2"][0][2].reshape(-1, 1)
        # damp2_* (K9)
        c["damp2_zeta1"] = np.array([z for z, _ in K.DAMP2_ZETA], float).reshape(1, -1)
        c["damp2_zeta2"] = np.array([z for _, z in K.DAMP2_ZETA], float).reshape(1, -1)
        c["damp2_f"] = .05
        c["damp2_E0"] = buf[nm]["damp2"][0][2].reshape(-1, 1)
        c["damp2_early"] = np.stack([p[0] for p in buf[nm]["damp2"]]).astype(np.float32)
        c["damp2_late"] = np.stack([p[1] for p in buf[nm]["damp2"]]).astype(np.float32)
        print(f"  {nm}: abs_regular={np.round(c['abs_regular'].ravel(),3)}", flush=True)

    savemat(out, dict(D=D, ECUT=K.ECUT, MU=K.MU, BETA2=K.BETA2,
                      fli_chaos_threshold=8., section_chaos_threshold=10.),
            do_compression=True, oned_as="column")
    print(f"[{mode}] salvo {out} ({out.stat().st_size/1e6:.2f} MB)")

    # round-trip: campos ORIGINAIS inalterados
    written = {"abs_f", "abs_maps", "abs_xs", "abs_vs", "abs_regular",
               "sec2_f", "sec2_pts", "sec2_fli", "sec2_E0",
               "damp2_zeta1", "damp2_zeta2", "damp2_f", "damp2_E0", "damp2_early", "damp2_late"}
    back = unwrap(loadmat(out, squeeze_me=False, struct_as_record=False)["D"])
    bad = 0
    for nm in cases:
        Csrc = unwrap(getattr(Dsrc, nm)); Cb = unwrap(getattr(back, nm))
        for f in Csrc._fieldnames:
            if f in written:      # campos 2.5 GL: adicionados/substituídos de propósito
                continue
            a = np.asarray(getattr(Csrc, f)); b = np.asarray(getattr(Cb, f))
            if a.shape != b.shape or not np.array_equal(a, b):
                print(f"  ❌ campo 1.5 GL alterado: {nm}.{f}  {a.shape}->{b.shape}"); bad += 1
    print("round-trip dos campos 1.5 GL:", "EXATO ✅" if bad == 0 else f"{bad} divergências ❌")


if __name__ == "__main__":
    main()
