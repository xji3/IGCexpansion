"""Write the results report (LaTeX) from report_data.json; compile with xelatex.
  python report/gen_report.py <report_data.json> <out.tex>"""
import json, math, sys

D = json.load(open(sys.argv[1]))
OUT = sys.argv[2]
cells = ['3', '10', '50', '100', '200', '300', '400', '500']
NAME = {'naive': 'Wald, naive', 'sandwich': 'Wald, sandwich', 'clr_global': r'CLR, global $\lambda$',
        'clr_local': r'CLR, local $\lambda$', 'pb': 'Bootstrap percentile'}


def P(v):
    """an endpoint on the psi scale, shown as 1/p"""
    if v is None:
        return r'$\infty$'
    x = math.exp(v)
    if x >= 1e5:
        e = int(math.log10(x))
        return f'${x / 10 ** e:.1f}{{\\times}}10^{{{e}}}$'
    if x >= 100:
        return f'{x:,.0f}'
    if x >= 10:
        return f'{x:.0f}'
    if x < 0.1:
        return '$<$0.1'
    return f'{x:.1f}'


def I(iv):
    return '--' if iv is None else f'[{P(iv[0])}, {P(iv[1])}]'


def B(s, bold):
    return r'\textbf{' + s + '}' if bold else s


def pval(p):
    if p is None:
        return '--'
    if p < 1e-4:
        return r'$<10^{-4}$'
    if p < 1e-3:
        return f'{p:.1e}'
    return f'{p:.3f}' if p < 0.1 else f'{p:.2f}'


L = []
MISS = r' $\times$'
w = L.append

# ---- point estimates
pe = []
for c in cells:
    v = D['cells'][c]
    pe.append(f"{c} & {v['n']} & {v['psi_hat_median_invp']:.1f} & [{v['psi_hat_q25_invp']:.1f}, {v['psi_hat_q75_invp']:.1f}] & "
              f"{v['frac_psi_hat_below']:.2f} & {v['boundary_reject']:.2f} \\\\")

# ---- main table: coverage / mean width / mean IS / median IS
meths = ['naive', 'sandwich', 'clr_global', 'clr_local']
main = []
for c in cells:
    ms = D['cells'][c]['methods']
    best_mean = min(ms[m]['is_mean'] for m in meths)
    best_med = min(ms[m]['is_med'] for m in meths)
    row = [c]
    for m in meths:
        q = ms[m]
        row.append(f"{q['cover']:.2f} & {q['width_mean']:.2f} & {B('%.2f' % q['is_mean'], q['is_mean'] == best_mean)} & "
                   f"{B('%.2f' % q['is_med'], q['is_med'] == best_med)}")
    main.append(' & '.join(row) + r' \\')

# ---- sandwich vs local split
split = []
for c in cells:
    ms = D['cells'][c]['methods']
    a, b = ms['sandwich'], ms['clr_local']
    split.append(f"{c} & {a['cover']:.2f} & {a['width_mean']:.2f} & {a['penalty_mean']:.2f} & "
                 f"{b['cover']:.2f} & {b['width_mean']:.2f} & {b['penalty_mean']:.2f} \\\\")

# ---- bootstrap pooled
pm = meths + ['pb']
pool = D['pb_pooled']
best_mean = min(pool[m]['is_mean'] for m in pm)
best_med = min(pool[m]['is_med'] for m in pm)
pooled = []
for m in pm:
    q = pool[m]
    pooled.append(f"{NAME[m]} & {round(q['cover'] * q['n'])}/{q['n']} & {q['width_mean']:.2f} & {q['width_med']:.2f} & "
                  f"{B('%.2f' % q['is_mean'], q['is_mean'] == best_mean)} & {B('%.2f' % q['is_med'], q['is_med'] == best_med)} & "
                  f"{q['miss_low']} / {q['miss_high']} \\\\")

# ---- bootstrap per case
cases = []
for e in D['pb_cases']:
    x = math.log(e['cell'])
    pb = e['pb']
    miss = not (pb[0] <= x <= (pb[1] if pb[1] is not None else 1e9))
    cases.append(f"{e['cell']} & {e['rep']} & {P(e['psi_hat'])} & {I(e['sandwich'])} & {I(e['clr_global'])} & "
                 f"{I(e['clr_local'])} & {I(e['pb'])}{MISS if miss else ''} \\\\")

# ---- real data
rd = []
for r in D['realdata']:
    pb = r['pb']
    pbs = '--' if pb is None else f"[{P(pb['lo'])}, {P(pb['hi'])}]"
    med = '--' if pb is None or pb.get('median') is None else P(pb['median'])
    lb = '--' if pb is None else str(pb.get('n_at_lower_bound'))
    name = r['name'].replace('_', r'\_') + (r'$^{a}$' if r['name'] == 'Harpak_introns' else '')
    rd.append(f"{name} & {P(r['psi_hat'])} & {pval(r['boundary_p'])} & {I(r['sandwich'])} & {I(r['clr_local'])} & "
              f"{pbs} & {med} & {lb} \\\\")

import os
T = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PSIGC_RESULTS_body.tex')).read()
for k, v in {'%%PE%%': pe, '%%MAIN%%': main, '%%SPLIT%%': split, '%%POOLED%%': pooled, '%%CASES%%': cases,
             '%%REAL%%': rd}.items():
    T = T.replace(k, '\n'.join(v))
open(OUT, 'w').write(T)
print('written', OUT)
