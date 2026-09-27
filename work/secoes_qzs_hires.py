"""Seções de Poincaré do QZS-ADV recalculadas e plotadas 100% em Python.

Porte fiel do MATLAB (qzs_integrar.m + qzs_config_beta_286.m + plot_analises_beta_286.m):
mesmo sistema, mesmo integrador (Yoshida 4a ordem simplético + mapa tangente para o FLI),
mesma classificação regular/caótico (ΔFLI = FLI(2t)-FLI(t) ≤ τ) e mesmos parâmetros.

Diferença proposital: as órbitas são plotadas com pontos grandes e coloridos por energia,
em alta resolução, corrigindo o problema do MATLAB (pontos finos e pequenos, difíceis de ver).

Validação: os percentuais ΔFLI>1 reproduzem o MATLAB órbita a órbita
  monoestável 16/31/12/42 %, poços rasos 17/25/47/39 %, poços profundos 2/3/9/12 %
  (com absorvedor, β=0.35; f = 0.002/0.01/0.05/0.15).

Uso:
  python secoes_qzs_hires.py                      # com absorvedor, β=0.35 (K8)
  python secoes_qzs_hires.py --sem-absorvedor     # sem absorvedor (K2)
  python secoes_qzs_hires.py --beta 0.25          # outra sintonia
"""
import argparse
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------ configuração
OMEGA = 0.35
MU = 0.1
TAU = 1.0
NORB = 80
PSECOES = 1500
NS_SECOES = round(200 / OMEGA)            # 571 passos por período de excitação
FSECOES = [0.0, 0.01, 0.05, 0.15]   # 1a coluna: f=0 (nao perturbado, referencia KAM)

CASES = {
    'monostable':    dict(eta=2/3.1, xw=1.4, vw=1.4, x0max=1.25, rotulo='Monostable (QZS)'),
    'shallow_wells': dict(eta=0.60,  xw=1.4, vw=1.4, x0max=1.30, rotulo='Shallow wells'),
    'deep_wells':    dict(eta=0.30,  xw=1.9, vw=1.5, x0max=1.75, rotulo='Deep wells'),
}
ORDEM = ['monostable', 'shallow_wells', 'deep_wells']


def U(x, eta):
    return 1.55 * x**2 - 3 * np.sqrt((1.5 * eta)**2 + x**2)


def xmin(eta):
    return 1.5 * np.sqrt(max((2 / 3.1)**2 - eta**2, 0.0))


# ------------------------------------------------------------------ integrador
def integrar(x0, v0, eta, beta, f, modelo, periodos, ns, snaps):
    """Yoshida 4 + mapa tangente (ζ=0). Retorna pts (periodos×2×n), growth e regular.

    Reproduz o ramo conservativo de qzs_integrar.m.
    """
    x = np.asarray(x0, float).copy()
    v = np.asarray(v0, float).copy()
    n = x.size
    x2 = x.copy(); v2 = v.copy()
    a2 = (1.5 * eta)**2
    b2 = beta**2; mb = MU * b2
    if modelo == 1:
        b2 = 0.0; mb = 0.0
    T = 2 * np.pi / OMEGA
    h = T / ns
    w1 = 1 / (2 - 2**(1/3)); w0 = -2**(1/3) / (2 - 2**(1/3))
    cy = np.array([w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]) * h
    dy = np.array([w1, w0, w1]) * h

    # força pré-calculada (t acumula em cy, como no MATLAB)
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
    snaps = list(snaps)
    fli = np.zeros((len(snaps), n))
    pts = np.zeros((periodos, 2, n), np.float32)

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
        pts[c-1, 0, :] = x.astype(np.float32)
        pts[c-1, 1, :] = v.astype(np.float32)

    growth = fli[1, :] - fli[0, :]
    return pts, growth, growth <= TAU


def secao_ci(nm):
    c = CASES[nm]
    pos = np.linspace(0.03, c['x0max'], NORB // 2)
    xline = np.concatenate([-pos[::-1], pos])
    return xline, np.zeros_like(xline)


def calcular(beta, modelo):
    dados = {}
    for nm in ORDEM:
        eta = CASES[nm]['eta']
        xline, vline = secao_ci(nm)
        E0 = U(xline, eta) - U(xmin(eta), eta)
        for f in FSECOES:
            t0 = time.time()
            pts, growth, regular = integrar(xline, vline, eta, beta, f, modelo,
                                            PSECOES, NS_SECOES, [PSECOES//2, PSECOES])
            dados[(nm, f)] = dict(pts=pts, regular=regular, E0=E0)
            print('  %-14s f=%-6g  ΔFLI>1=%2.0f%%  (%.0fs)'
                  % (nm, f, 100*np.mean(~regular), time.time()-t0), flush=True)
    return dados


# ------------------------------------------------------------------ figura
def _cmap_energia():
    rgb = np.array([[43,26,111],[59,91,219],[28,159,214],[32,178,170],
                    [108,194,74],[224,195,0],[240,140,0],[215,38,61]]) / 255
    return LinearSegmentedColormap.from_list('energia', rgb, N=512)


def plotar(dados, beta, modelo, base, dpi):
    cmap_e = _cmap_energia()
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
    fig, axs = plt.subplots(3, 4, figsize=(15, 11.5), dpi=dpi)
    fig.patch.set_facecolor('white')
    sc = None
    for r, nm in enumerate(ORDEM):
        c0 = CASES[nm]; eta = c0['eta']; xw, vw = c0['xw'], c0['vw']
        Um = U(xmin(eta), eta); depth = U(0., eta) - Um
        xg = np.linspace(-xw, xw, 400); vg = np.linspace(-vw, vw, 400)
        XG, VG = np.meshgrid(xg, vg); HG = 0.5*VG**2 + U(XG, eta) - Um
        for c, f in enumerate(FSECOES):
            ax = axs[r, c]; ax.set_facecolor('white')
            d = dados[(nm, f)]; pts = d['pts']; reg = d['regular']; E0 = d['E0']
            chaos = ~reg
            if chaos.any():
                ax.scatter(pts[:, 0, chaos].ravel(), pts[:, 1, chaos].ravel(),
                           s=1.4, c='#8a8a8a', linewidths=0, alpha=0.35, rasterized=True)
            if reg.any():
                idx = np.where(reg)[0]
                col = np.repeat(np.clip(E0[idx], 0, 0.8)[None, :], pts.shape[0], axis=0)
                sc = ax.scatter(pts[:, 0, idx].ravel(), pts[:, 1, idx].ravel(),
                                s=2.2, c=col.ravel(), cmap=cmap_e, vmin=0, vmax=0.8,
                                linewidths=0, alpha=0.75, rasterized=True)
            if depth > 1e-9:
                ax.contour(XG, VG, HG, [depth], colors='0.15', linewidths=0.6, linestyles='--')
            ax.set_xlim(-xw, xw); ax.set_ylim(-vw, vw); ax.set_aspect('equal', 'box')
            ax.tick_params(labelsize=13)
            ax.text(0.04, 0.96, r'$\Delta$FLI>1: %.0f%%' % (100*np.mean(chaos)),
                    transform=ax.transAxes, va='top', ha='left', fontsize=14,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.5', alpha=0.9))
            if r == 0:
                ax.set_title('f = 0 (KAM)' if f == 0 else ('f = %g' % f),
                             fontsize=19, fontweight='bold')
            if r == 2:
                ax.set_xlabel('X', fontsize=17)
            else:
                ax.set_xticklabels([])
            if c == 0:
                ax.set_ylabel('V = dX/dt', fontsize=17)
            else:
                ax.set_yticklabels([])
        axs[r, 0].annotate('%s ($\\eta$=%.4g)' % (c0['rotulo'], eta),
                           xy=(-0.47, 0.5), xycoords='axes fraction', rotation=90,
                           ha='center', va='center', fontsize=16)

    if modelo == 2:
        tit = r'Poincaré sections (X, V) — with absorber;  $\Omega$=%.2f,  $\beta$=%.2f' % (OMEGA, beta)
    else:
        tit = r'Poincaré sections (X, V) — without absorber;  $\Omega$=%.2f' % OMEGA
    fig.suptitle(tit, fontsize=19, y=0.985)
    cbax = fig.add_axes([0.30, 0.055, 0.40, 0.018])
    cb = fig.colorbar(sc, cax=cbax, orientation='horizontal')
    cb.set_label(r'initial energy of the primary  $E_0 - U_{\min}$  (regular orbits)', fontsize=15)
    cb.ax.tick_params(labelsize=13)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.94, bottom=0.135, wspace=0.08, hspace=0.10)
    fig.savefig(base + '.png', dpi=dpi, facecolor='white')
    fig.savefig(base + '.pdf', facecolor='white')
    print('figuras salvas:', base + '.png', '/', base + '.pdf')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--beta', type=float, default=0.35)
    ap.add_argument('--sem-absorvedor', action='store_true', help='modelo 1.5 GL (K2)')
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    modelo = 1 if a.sem_absorvedor else 2
    tag = 'sem_absorvedor' if modelo == 1 else ('com_absorvedor_beta%.2f' % a.beta)
    base = a.out or ('secoes_%s_hires' % tag)
    print('Calculando seções (%s)...' % tag)
    dados = calcular(a.beta, modelo)
    plotar(dados, a.beta, modelo, base, a.dpi)


if __name__ == '__main__':
    main()
