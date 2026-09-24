"""Rendering/analysis helpers for the KAM study (work/kam_compute.py).

Currently exposes `resonances`, used by export_kam_matlab.py to locate the
p:q resonances of the unperturbed frequency map omega(E) of the reduced
1.5-DOF model. A resonance is a point where omega(E)/Omega = p/q with the
driving frequency Omega = 1, i.e. omega(E) = p/q.

omega(E) is sampled on two physical branches, tagged by `over` in the
omega_<case>.npz files produced by kam_compute.omega_curve:
  over = 0  motion inside the wells (present only when depth > 0)
  over = 1  motion over the barrier / single well
On each branch omega(E) is monotonic (d omega/dE < 0 inside the wells,
d omega/dE > 0 above), so each reachable ratio p/q is crossed at most once
per branch; both branches can host the same ratio (e.g. p/q = 1/2 inside a
well and again over the barrier). Crossings are located by linear
interpolation of the sampled omega(E) curve.
"""
from math import gcd

import numpy as np


def resonances(d, qmax=4, Omega=1.0):
    """Return the p:q resonances of omega(E), one dict per crossing.

    Args:
        d: mapping with arrays ``dE`` (E - Umin), ``omega_quad`` (omega by
           quadrature) and ``over`` (branch tag), as stored in
           omega_<case>.npz.
        qmax: largest denominator q to consider (callers filter further;
              export keeps q <= 3).
        Omega: driving frequency (1 by construction of the model).

    Returns:
        List of dicts with keys ``p``, ``q``, ``rho`` (= p/q), ``dE`` (the
        interpolated energy gap of the crossing) and ``over`` (branch tag),
        sorted by branch then by dE.
    """
    dE = np.asarray(d["dE"], float).ravel()
    omega = np.asarray(d["omega_quad"], float).ravel()
    over = np.asarray(d["over"], float).ravel()

    out = []
    for branch in sorted(np.unique(over)):
        m = over == branch
        if m.sum() < 2:
            continue
        de_b = dE[m]
        w_b = omega[m]
        order = np.argsort(de_b)
        de_b = de_b[order]
        w_b = w_b[order]
        finite = np.isfinite(w_b)
        de_b = de_b[finite]
        w_b = w_b[finite]
        if de_b.size < 2:
            continue
        wlo, whi = float(np.min(w_b)), float(np.max(w_b))
        for q in range(1, qmax + 1):
            for p in range(1, int(np.floor(q * whi / Omega)) + 1):
                if gcd(p, q) != 1:
                    continue
                rho = p / q
                target = rho * Omega
                if target < wlo or target > whi:
                    continue
                # linear-interpolate every sign change of (omega - target)
                g = w_b - target
                sign = np.sign(g)
                for i in np.where(sign[:-1] * sign[1:] < 0)[0]:
                    x0, x1 = de_b[i], de_b[i + 1]
                    y0, y1 = g[i], g[i + 1]
                    de_c = x0 + (x1 - x0) * (-y0) / (y1 - y0)
                    out.append(dict(p=p, q=q, rho=rho, dE=float(de_c),
                                    over=float(branch)))
                # exact hits on a sample point
                for i in np.where(g == 0.0)[0]:
                    out.append(dict(p=p, q=q, rho=rho, dE=float(de_b[i]),
                                    over=float(branch)))
    out.sort(key=lambda r: (r["over"], r["dE"]))
    return out
