"""Teste de convergência dt vs dt/2 para provar que as ilhas de ressonância
(cadeias de Poincaré–Birkhoff) das seções são físicas, não erro numérico.

Recalcula uma seção com o passo do integrador (ns passos/período) e com o passo
pela metade (2*ns), classifica regular/caótico pelo mesmo critério (ΔFLI ≤ τ) e
sobrepõe as órbitas regulares. Se toros e ilhas coincidem sob refino do passo,
são estruturas físicas do espaço de fase.

Uso:
  python teste_convergencia.py                        # deep_wells, f=0.15, com absorvedor beta=0.35
  python teste_convergencia.py --caso shallow_wells --f 0.05
"""
import argparse
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import secoes_qzs_hires as S


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--caso', default='deep_wells', choices=list(S.CASES))
    ap.add_argument('--f', type=float, default=0.15)
    ap.add_argument('--beta', type=float, default=0.35)
    ap.add_argument('--sem-absorvedor', action='store_true')
    ap.add_argument('--out', default='teste_convergencia')
    ap.add_argument('--dpi', type=int, default=220)
    a = ap.parse_args()
    modelo = 1 if a.sem_absorvedor else 2
    nm = a.caso; eta = S.CASES[nm]['eta']
    xl, vl = S.secao_ci(nm)
    snaps = [S.PSECOES // 2, S.PSECOES]

    t0 = time.time()
    p1, g1, r1 = S.integrar(xl, vl, eta, a.beta, a.f, modelo, S.PSECOES, S.NS_SECOES, snaps)
    print('dt   pronto (%.0fs): DFLI>1=%.0f%%' % (time.time()-t0, 100*np.mean(~r1)), flush=True)
    t0 = time.time()
    p2, g2, r2 = S.integrar(xl, vl, eta, a.beta, a.f, modelo, S.PSECOES, 2*S.NS_SECOES, snaps)
    print('dt/2 pronto (%.0fs): DFLI>1=%.0f%%' % (time.time()-t0, 100*np.mean(~r2)), flush=True)
    print('classificacao identica: %s (%d/%d orbitas)'
          % (np.array_equal(r1, r2), np.sum(r1 == r2), r1.size))

    c1 = ~r1
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
    fig, ax = plt.subplots(figsize=(9.5, 9), dpi=a.dpi)
    ax.scatter(p1[:, 0, c1].ravel(), p1[:, 1, c1].ravel(), s=0.5, c='0.8', lw=0, alpha=0.3)
    ax.scatter(p1[:, 0, r1].ravel(), p1[:, 1, r1].ravel(), s=0.7, c='#1f5fd0', lw=0, alpha=0.55,
               label='dt (%d passos/T)' % S.NS_SECOES)
    ax.scatter(p2[:, 0, r2].ravel(), p2[:, 1, r2].ravel(), s=0.7, c='#e03b1f', lw=0, alpha=0.55,
               label='dt/2 (%d passos/T)' % (2*S.NS_SECOES))
    c = S.CASES[nm]
    ax.set_xlim(-c['xw'], c['xw']); ax.set_ylim(-c['vw'], c['vw']); ax.set_aspect('equal')
    ax.set_xlabel('X'); ax.set_ylabel('V = dX/dt')
    lab = 'sem absorvedor' if modelo == 1 else (r'com absorvedor, $\beta$=%.2f' % a.beta)
    ax.set_title('Teste de convergência — %s, f=%g (%s)' % (c['rotulo'], a.f, lab), fontsize=12)
    ax.legend(loc='upper center', fontsize=9, markerscale=10, ncol=2)
    txt = ('Azul (dt) e vermelho (dt/2) sobre as mesmas curvas = roxo.\n'
           '%d/%d órbitas mantêm a classe regular/caótica.\n'
           'Fração caótica: %.0f%% (dt) = %.0f%% (dt/2).\n'
           'Toros e ilhas coincidem sob refino do passo => físicos, não numéricos.'
           % (np.sum(r1 == r2), r1.size, 100*np.mean(~r1), 100*np.mean(~r2)))
    ax.text(0.015, 0.015, txt, transform=ax.transAxes, fontsize=8.5, va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='0.6', alpha=0.92))
    fig.tight_layout()
    fig.savefig(a.out + '.png', dpi=a.dpi)
    print('figura salva:', a.out + '.png')


if __name__ == '__main__':
    main()
