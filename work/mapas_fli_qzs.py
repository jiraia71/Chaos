"""Mapas de FLI do QZS-ADV (comparação K5) recalculados e plotados em Python.

Mesmo integrador simplético + mapa tangente de secoes_qzs_hires.py, mas sobre uma
grade 200×200 de condições iniciais (X0, V0). Cor = FLI(400T); percentual "regular*"
= fração da área com ΔFLI = FLI(400T)-FLI(200T) ≤ τ, dentro de H0 - Umin ≤ ecut.

K5: poços rasos, linha de cima "sem absorvedor" (1.5 GL), linha de baixo
"com absorvedor" (2.5 GL). Colunas f = 0.01/0.05/0.15.
"""
import argparse
import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OMEGA = 0.35
MU = 0.1
TAU = 1.0
ECUT = 0.5
NX = 200
PMAPAS = 400
NS_MAPAS = round(100 / OMEGA)     # 286
FMAPAS = [0.01, 0.05, 0.15]

# K5 é definida para poços rasos
CASO = dict(nome='shallow_wells', eta=0.60, xw=1.4, vw=1.4, rotulo='Poços rasos')


def U(x, eta):
    return 1.55 * x**2 - 3 * np.sqrt((1.5 * eta)**2 + x**2)


def xmin(eta):
    return 1.5 * np.sqrt(max((2 / 3.1)**2 - eta**2, 0.0))


def mapa(eta, beta, f, modelo, nx=NX, periodos=PMAPAS, ns=NS_MAPAS):
    """FLI(400T) e growth sobre a grade; retorna (xs, vs, fli_map, growth_map, inside)."""
    xs = np.linspace(-CASO['xw'], CASO['xw'], nx)
    vs = np.linspace(-CASO['vw'], CASO['vw'], nx)
    X, V = np.meshgrid(xs, vs)
    x = X.ravel().copy(); v = V.ravel().copy()
    n = x.size
    x2 = x.copy(); v2 = v.copy()
    a2 = (1.5 * eta)**2
    b2 = beta**2; mb = MU * b2
    if modelo == 1:
        b2 = 0.0; mb = 0.0
    T = 2 * np.pi / OMEGA; h = T / ns
    w1 = 1 / (2 - 2**(1/3)); w0 = -2**(1/3) / (2 - 2**(1/3))
    cy = np.array([w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]) * h
    dy = np.array([w1, w0, w1]) * h
    force = np.zeros((ns, 3)); t = 0.0
    for s in range(ns):
        for k in range(4):
            t += cy[k]
            if k < 3:
                force[s, k] = f * np.cos(OMEGA * t)
    if modelo == 1:
        dx = np.full(n, 1/np.sqrt(2)); dv = dx.copy(); dx2 = np.zeros(n); dv2 = np.zeros(n)
    else:
        dx = np.full(n, 0.5); dv = dx.copy(); dx2 = dx.copy(); dv2 = dx.copy()
    logs = np.zeros(n)
    snaps = [periodos // 2, periodos]
    fli = np.zeros((2, n))
    for c in range(1, periodos + 1):
        for s in range(ns):
            for k in range(4):
                x += cy[k] * v
                if modelo == 2:
                    x2 += cy[k] * v2
                dx += cy[k] * dv
                if modelo == 2:
                    dx2 += cy[k] * dv2
                if k == 3:
                    break
                rad = np.sqrt(a2 + x**2)
                du = 3.1 * x - 3 * x / rad
                v += dy[k] * (force[s, k] - du - mb * (x - x2))
                if modelo == 2:
                    v2 += dy[k] * b2 * (x - x2)
                curv = 3.1 - 3 * a2 / rad**3
                dv += dy[k] * (-(curv + mb) * dx + mb * dx2)
                if modelo == 2:
                    dv2 += dy[k] * b2 * (dx - dx2)
        normw = np.sqrt(dx**2 + dv**2 + dx2**2 + dv2**2)
        logs += np.log10(normw)
        dx /= normw; dv /= normw; dx2 /= normw; dv2 /= normw
        if c in snaps:
            fli[snaps.index(c), :] = logs
    growth = (fli[1] - fli[0]).reshape(nx, nx)
    Um = U(xmin(eta), eta)
    H0 = 0.5 * V**2 + U(X, eta) - Um
    inside = H0 <= ECUT
    return xs, vs, fli[1].reshape(nx, nx), growth, inside


def _cmap_fli():
    # azul -> vermelho (paleta FLI de referência)
    rgb = np.array([[8,29,88],[34,94,168],[29,145,192],[65,182,196],[127,205,187],
                    [199,233,180],[255,255,204],[254,204,92],[253,141,60],
                    [240,59,32],[189,0,38]]) / 255
    return LinearSegmentedColormap.from_list('fli', rgb, N=512)


CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'dados', 'cache_mapas_fli')


def get_map(eta, beta, f, modelo, cache=CACHE):
    """Calcula (ou carrega do cache) um mapa. modelo=1 nao depende de beta."""
    os.makedirs(cache, exist_ok=True)
    bkey = 0.0 if modelo == 1 else beta      # sem absorvedor independe de beta
    path = os.path.join(cache, 'map_m%d_b%.3g_f%.4g.npz' % (modelo, bkey, f))
    if os.path.exists(path):
        d = np.load(path)
        return d['xs'], d['vs'], d['M'], d['G'], d['inside']
    t0 = time.time()
    xs, vs, M, G, inside = mapa(eta, beta, f, modelo)
    np.savez_compressed(path, xs=xs, vs=vs, M=M.astype(np.float32),
                        G=G.astype(np.float32), inside=inside)
    print('  computed m%d f=%-5g (%.0fs)' % (modelo, f, time.time()-t0), flush=True)
    return xs, vs, M, G, inside


def gerar(beta, base, dpi, cache=CACHE):
    cmap = _cmap_fli()
    eta = CASO['eta']; Um = U(xmin(eta), eta); depth = U(0., eta) - Um
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
    fig, axs = plt.subplots(2, 3, figsize=(12, 9.0), dpi=dpi)
    fig.patch.set_facecolor('white')
    im = None
    for r, modelo in enumerate([1, 2]):
        for c, f in enumerate(FMAPAS):
            xs, vs, M, G, inside = get_map(eta, beta, f, modelo, cache)
            reg = np.mean(G[inside] <= TAU)
            ax = axs[r, c]
            im = ax.imshow(np.clip(M, 0, 20), origin='lower',
                           extent=[xs[0], xs[-1], vs[0], vs[-1]], aspect='equal',
                           cmap=cmap, vmin=0, vmax=20)
            X, V = np.meshgrid(xs, vs); H = 0.5*V**2 + U(X, eta) - Um
            ax.contour(X, V, H, [ECUT], colors='w', linewidths=0.6, linestyles=':')
            if depth > 1e-9:
                ax.contour(X, V, H, [depth], colors='w', linewidths=0.5, linestyles='--')
            ax.text(0.04, 0.96, 'regular* %.1f%%' % (100*reg), transform=ax.transAxes,
                    va='top', ha='left', color='w', fontsize=14,
                    bbox=dict(boxstyle='round,pad=0.25', fc=(.08, .15, .25), ec='none', alpha=0.85))
            if r == 0:
                ax.set_title('f = %g' % f, fontsize=18, fontweight='bold')
            if r == 1:
                ax.set_xlabel(r'$X_0$', fontsize=17)
            else:
                ax.set_xticklabels([])
            if c == 0:
                ax.set_ylabel(r'$V_0$', fontsize=17)
            else:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=13)
        lab = 'Without absorber' if modelo == 1 else (r'With absorber ($\beta$=%.2f)' % beta)
        axs[r, 0].annotate(lab, xy=(-0.36, 0.5), xycoords='axes fraction', rotation=90,
                           ha='center', va='center', fontsize=16)
    fig.suptitle(r'FLI maps — shallow wells;  $\Omega$=%.2f,  $\beta$=%.2f' % (OMEGA, beta),
                 fontsize=18, y=0.975)
    cbax = fig.add_axes([0.30, 0.085, 0.40, 0.020])
    cb = fig.colorbar(im, cax=cbax, orientation='horizontal')
    cb.set_label(r'FLI(%dT), $\log_{10}$ (scale capped at 20)' % PMAPAS, fontsize=15)
    cb.ax.tick_params(labelsize=13)
    fig.subplots_adjust(left=0.17, right=0.985, top=0.905, bottom=0.185, wspace=0.06, hspace=0.08)
    fig.savefig(base + '.png', dpi=dpi, facecolor='white')
    fig.savefig(base + '.pdf', facecolor='white')
    print('saved:', base + '.png', '/', base + '.pdf')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--beta', type=float, default=0.35)
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    base = a.out or ('fli_maps_beta%.2f' % a.beta)
    gerar(a.beta, base, a.dpi)


if __name__ == '__main__':
    main()
