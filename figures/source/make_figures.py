"""Rebuild every manuscript panel from frozen closeout CSVs; no model fitting.

Run from any directory: python path/to/figures/source/make_figures.py
PDF is the Overleaf master; SVG is editable; PNG is a 400-dpi preview.
"""
from pathlib import Path
import json, hashlib, sys, importlib.metadata
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / 'figures' / 'inputs'
OUT = ROOT / 'paper' / 'figures'
PAIRS = ['ddec6_cm5', 'ddec6_bader', 'cm5_bader']
NAMES = ['DDEC6 − CM5', 'DDEC6 − Bader', 'CM5 − Bader']
BLUE, TEAL, ORANGE, GRAY = '#31688E', '#16857B', '#D47C38', '#7B8490'
SPLITS = [('original', BLUE, 'o', 'Structure folds'),
          ('blocked', TEAL, 's', 'Identifier/composition blocks')]
LEFT = ['C:O-attached', 'C:N-attached', 'H:O-bound', 'H:N-bound',
        'O:metal-adjacent', 'N:metal-adjacent']
LABELS = ['C–O vs C/H carbon', 'C–N vs C/H carbon', 'O–H vs C–H',
          'N–H vs C–H', 'O: metal vs no metal', 'N: metal vs no metal']
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10,
    'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .65, 'lines.linewidth': 1.2, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white',
    'axes.titlepad': 9, 'figure.constrained_layout.use': True,
    'figure.constrained_layout.w_pad': .07, 'figure.constrained_layout.h_pad': .07,
})
USED = {}

def read(name):
    USED[name] = True
    return pd.read_csv(IN / name)

def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{name}.{ext}', dpi=400,
                    metadata={'Creator': 'ChaMM frozen-table figure pipeline'} if ext == 'pdf' else None)
    plt.close(fig)

def panel(ax, letter):
    ax.text(0, 1.10, letter, transform=ax.transAxes, weight='bold',
            fontsize=11, ha='left', va='bottom')

def point(ax, value, lo, hi, y, color, marker='o', scale=1):
    if not np.isfinite([value, lo, hi]).all():
        return
    ax.errorbar(value * scale, y,
                xerr=[[max(0, value-lo)*scale], [max(0, hi-value)*scale]],
                fmt=marker, ms=3.7, color=color, ecolor=color,
                capsize=2, elinewidth=.9, markeredgewidth=.75, zorder=4)

def split_legend(ax, loc='best'):
    handles = [Line2D([], [], color=c, marker=m, linestyle='none', markersize=4,
                      label=label) for _, c, m, label in SPLITS]
    ax.figure.legend(handles=handles, frameon=False, loc='outside lower center', ncol=2)

def fig1():
    d = read('element_differences.csv')
    els = ['H', 'C', 'N', 'O', 'Mg', 'Al', 'Ca', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Zr', 'Cd', 'Pb']
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 4.9), sharey=True)
    counts = d.drop_duplicates('element').set_index('element').n_structures
    for j, (ax, pair) in enumerate(zip(axs, PAIRS)):
        g = d[d.pair.eq(pair)].set_index('element').reindex(els)
        for i, (_, r) in enumerate(g.iterrows()):
            color = ORANGE if bool(r.is_metal) else BLUE
            ax.plot([0, r.structure_balanced_mean_signed], [i, i], color=color, alpha=.55, lw=1.7)
            ax.scatter(r.structure_balanced_mean_signed, i, s=18, color=color, zorder=3)
        ax.axvline(0, color=GRAY, lw=.65)
        ax.set(title=NAMES[j], xlabel='Mean signed difference (e)')
        ax.grid(axis='x', color='#E6E9EC', lw=.6)
        panel(ax, chr(97+j))
    axs[0].set_yticks(range(len(els)), [f'{e}   {int(counts[e]):,}' for e in els])
    axs[0].invert_yaxis()
    axs[0].set_ylabel('Element and structure count')
    fig.legend(handles=[Line2D([], [], marker='o', linestyle='none', color=BLUE, label='Nonmetal'),
                           Line2D([], [], marker='o', linestyle='none', color=ORANGE, label='Metal')],
                   frameon=False, loc='outside lower center', ncol=2)
    save(fig, 'Figure_01_charge_landscape')

def fig2():
    d, inc = read('conditional_metrics.csv'), read('conditional_gains.csv')
    models = ['piecewise', 'charge_only', 'conditional_first', 'conditional_twohop']
    fig, axs = plt.subplots(2, 3, figsize=(7.3, 6.2), gridspec_kw={'height_ratios': [1.4, 1]})
    for j, pair in enumerate(PAIRS):
        ax = axs[0, j]
        for k, (split, color, marker, _) in enumerate(SPLITS):
            g = d[d.pair.eq(pair) & d.split.eq(split) & d.stratum.eq('all')].set_index('model').loc[models]
            x = np.arange(4)+(k-.5)*.11
            ax.errorbar(x, g.mae_structure, yerr=[g.mae_structure-g.ci_low, g.ci_high-g.mae_structure],
                        fmt=marker+'-', color=color, ms=3.5, capsize=2, lw=1)
        ax.set(title=NAMES[j], xticks=range(4),
               xticklabels=['12-bin\nmapping', 'Charge\nonly', 'First\nshell', 'Two-hop\npaths'])
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', color='#E6E9EC', lw=.6)
        if j == 0: ax.set_ylabel('Structure-balanced MAE (e)')
        panel(ax, chr(97+j))
        bx = axs[1, j]
        for k, (split, color, marker, _) in enumerate(SPLITS):
            for y, pop in enumerate(['all', 'has_mofkey', 'nonsingleton']):
                r = inc[inc.pair.eq(pair) & inc.split.eq(split) & inc.stratum.eq(pop)
                        & inc.comparison.eq('structure_given_q_A')].iloc[0]
                point(bx, r.improvement_structure, r.ci_low, r.ci_high, y+(k-.5)*.19, color, marker, 1000)
        bx.axvline(0, color=GRAY, lw=.7)
        bx.set_yticks(range(3), ['All: 9,171', 'MOFkey: 5,408', 'Nonsingleton: 1,414'] if j == 0 else ['', '', ''])
        bx.invert_yaxis()
        bx.set_xlabel('Gain over charge only (10⁻³ e)')
        panel(bx, chr(100+j))
    split_legend(axs[1, 2], 'lower right')
    save(fig, 'Figure_02_conditional_evidence')

def contrast_grid(d, name, splits=True):
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 6.8), sharey=True)
    for row, quantity in enumerate(['calibrated_signed', 'calibrated_absolute']):
        for j, pair in enumerate(PAIRS):
            ax = axs[row, j]
            for k, (split, color, marker, _) in enumerate(SPLITS):
                for y, left in enumerate(LEFT):
                    r = d[d.pair.eq(pair) & d.split.eq(split) & d.left.eq(left) & d.quantity.eq(quantity)].iloc[0]
                    point(ax, r.effect, r.ci_low, r.ci_high, y+(k-.5)*.18, color, marker)
            ax.axvline(0, color=GRAY, lw=.7)
            ax.set_title(NAMES[j])
            ax.grid(axis='x', color='#E6E9EC', lw=.6)
            ax.set_xlabel(('Signed' if row == 0 else 'Absolute')+' departure contrast (e)')
            if j == 0:
                ns = d[d.pair.eq(pair) & d.split.eq('blocked') & d.quantity.eq(quantity)].set_index('left').n_structures
                ax.set_yticks(range(6), [f'{label}\n(n = {int(ns[left]):,})' for label, left in zip(LABELS, LEFT)])
            panel(ax, chr(97+row*3+j))
    axs[0, 0].invert_yaxis()
    split_legend(axs[1, 2], 'lower right')
    save(fig, name)

def fig4():
    d = read('conditional_gains.csv')
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.1), sharey=True)
    for j, pair in enumerate(PAIRS):
        ax = axs[j]
        for k, (split, color, marker, _) in enumerate(SPLITS):
            for y, pop in enumerate(['all', 'nonmetal', 'metal']):
                r = d[d.pair.eq(pair) & d.split.eq(split) & d.stratum.eq(pop)
                      & d.comparison.eq('structure_given_q_A')].iloc[0]
                point(ax, r.improvement_structure, r.ci_low, r.ci_high, y+(k-.5)*.18, color, marker, 1000)
        ax.axvline(0, color=GRAY, lw=.8)
        ax.set(title=NAMES[j], xlabel='Gain over charge only (10⁻³ e)')
        ax.grid(axis='x', color='#E6E9EC', lw=.6)
        panel(ax, chr(97+j))
    axs[0].set_yticks(range(3), ['All atoms', 'Nonmetal atoms', 'Metal atoms'])
    axs[0].invert_yaxis()
    split_legend(axs[-1], 'upper right')
    save(fig, 'Figure_04_population_limits')

def fig5():
    changes, edges = read('nonmetal_changes.csv'), read('nonmetal_edges.csv')
    g = changes.groupby('scale')[['first_shell_changed', 'twohop_changed', 'environment_changed']].sum()
    fig = plt.figure(figsize=(7.2, 4.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
    ax, bx, cx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    x = np.arange(3)
    for k, (col, color, label) in enumerate([('first_shell_changed', BLUE, 'First shell'),
                                            ('twohop_changed', TEAL, 'Two-hop paths')]):
        v = g[col].reindex([1.2, 1.25, 1.3], fill_value=0).to_numpy()
        ax.bar(x+(k-.5)*.24, v, width=.23, color=color, label=label)
        for xx, vv in zip(x+(k-.5)*.24, v): ax.text(xx, vv+.7, str(int(vv)), ha='center', fontsize=7)
    ax.set(xticks=x, xticklabels=['1.20', '1.25\naccepted', '1.30\nstress'], ylabel='Atoms with changed descriptors', ylim=(0, 48))
    ax.legend(frameon=False, loc='upper left'); panel(ax, 'a')
    # Geometry is reconstructed from actual recorded O-H lengths and O-H-O angle.
    e = edges[edges.accepted_H_cn.eq(0)].iloc[0]
    same = edges[edges.qmof_id.eq(e.qmof_id) & (edges.i.eq(e.i) | edges.j.eq(e.i))]
    distances = sorted(same.distance.tolist())[:2]
    angle = float(e.H_angle_degree)
    left = np.array([-distances[0], 0.])
    right = np.array([-distances[1]*np.cos(np.deg2rad(angle)), distances[1]*np.sin(np.deg2rad(angle))])
    for v in [left, right]: bx.plot([0,v[0]], [0,v[1]], '--', color=ORANGE, lw=1.6)
    bx.scatter([left[0],right[0]], [left[1],right[1]], s=140, color='#BA5661', zorder=3)
    bx.scatter(0, 0, s=55, color='#B8BDC3', edgecolor=GRAY, zorder=4)
    for v, label in [(left, 'O'), (right, 'O'), (np.zeros(2), 'H')]:
        bx.text(v[0], v[1], label, ha='center', va='center', fontsize=7, color='white' if label=='O' else '#20262C')
    for v, distance in zip([left,right], distances): bx.text(v[0]/2, .19, f'{distance:.4f} Å', ha='center', fontsize=7)
    bx.text(0,-.32,f'O–H–O: {angle:.3f}°',ha='center',fontsize=8)
    bx.text(0,-.66,'H degree: 0 → 2 at the wider cutoff',ha='center',fontsize=7)
    bx.set(xlim=(-1.7,1.7),ylim=(-.85,.55),title='A shared-proton ambiguity'); bx.set_aspect('equal'); bx.axis('off'); panel(bx,'b')
    # Fixed models from the earlier metal-only follow-up: deliberately separate from conditional tests.
    sens = read('metal_graph_sensitivity.csv')
    for k,(variant,color,marker,label) in enumerate([('shift',TEAL,'o','Metal taper +0.10 Å'),
                                                    ('cutoff',ORANGE,'s','Original radius-cutoff stress')]):
        for y,pair in enumerate(PAIRS):
            r=sens[sens.variant.eq(variant)&sens.pair.eq(pair)&sens.model.eq('metal_tree_second')].iloc[0]
            point(cx,r.alternative_minus_accepted,r.ci_low,r.ci_high,y+(k-.5)*.19,color,marker,1000)
    cx.axvline(0,color=GRAY,lw=.7)
    cx.set_yticks(range(3),NAMES);cx.invert_yaxis()
    cx.set_xlabel('Alternative − accepted MAE, earlier metal-only model (10⁻³ e)')
    fig.legend(handles=[Line2D([],[],color=TEAL,marker='o',linestyle='none',label='Metal taper +0.10 Å'),
                       Line2D([],[],color=ORANGE,marker='s',linestyle='none',label='Original radius-cutoff stress')],frameon=False,loc='outside lower center',ncol=2)
    cx.grid(axis='x',color='#E6E9EC',lw=.6);panel(cx,'c')
    save(fig,'Figure_05_connectivity_limits')

def supplement():
    d=read('calibration_curvature.csv')
    fig,axs=plt.subplots(2,3,figsize=(7.2,4.7),sharey='row')
    for i,el in enumerate(['C','O']):
        for j,pair in enumerate(PAIRS):
            ax=axs[i,j];g=d[d.element.eq(el)&d.pair.eq(pair)].sort_values('qA_mean')
            ax.plot(g.qA_mean,g.affine_signed_residual,'o-',color=ORANGE,ms=3,label='Affine')
            ax.plot(g.qA_mean,g.nonlinear_signed_residual,'s-',color=TEAL,ms=3,label='12-bin piecewise')
            ax.axhline(0,color=GRAY,lw=.6);ax.set(title=f'{el}: {NAMES[j]}',xlabel='Source-scheme charge (e)')
            if j==0:ax.set_ylabel('Mean signed residual (e)')
            panel(ax,chr(97+i*3+j))
    axs[0,2].legend(frameon=False);save(fig,'Figure_S01_calibration')
    d=read('conditional_metrics.csv')
    fig,axs=plt.subplots(1,3,figsize=(7.2,2.7))
    for j,pair in enumerate(PAIRS):
        ax=axs[j];g=d[d.split.eq('blocked')&d.pair.eq(pair)&d.stratum.eq('all')].set_index('model')
        for k,model in enumerate(['charge_only','conditional_first','conditional_twohop']):
            r=g.loc[model];ax.plot([r.mae_atom,r.mae_structure],[k,k],color=GRAY)
            ax.scatter(r.mae_atom,k,color=BLUE,s=18);ax.scatter(r.mae_structure,k,color=TEAL,s=18,marker='s')
        ax.set(title=NAMES[j],xlabel='Held-out MAE (e)',yticks=range(3),yticklabels=['Charge only','First shell','Two-hop paths'] if j==0 else ['','',''])
        ax.invert_yaxis();panel(ax,chr(97+j))
    fig.legend(loc='outside lower center',ncol=2,handles=[Line2D([],[],color=BLUE,marker='o',linestyle='none',label='Atom weighted'),Line2D([],[],color=TEAL,marker='s',linestyle='none',label='Structure balanced')],frameon=False)
    save(fig,'Figure_S02_estimands')
    d=read('fold_gains.csv')
    fig,axs=plt.subplots(1,3,figsize=(7.2,3.0),sharey=True)
    for j,pair in enumerate(PAIRS):
        ax=axs[j]
        for k,(split,color,marker,_) in enumerate(SPLITS):
            g=d[d.split.eq(split)&d.pair.eq(pair)&d.stratum.eq('all')&d.comparison.eq('structure_given_q_A')].sort_values('fold')
            ax.plot(g.improvement_structure*1000,np.arange(5)+(k-.5)*.15,marker+'-',color=color,ms=4)
        ax.axvline(0,color=GRAY,lw=.7);ax.set(title=NAMES[j],xlabel='Gain over charge only (10⁻³ e)',yticks=range(5),yticklabels=[str(i+1) for i in range(5)])
        panel(ax,chr(97+j))
    axs[0].set_ylabel('Held-out fold');axs[0].invert_yaxis();split_legend(axs[-1]);save(fig,'Figure_S03_folds')
    atlas=read('compact_atlas.csv')
    envs=['C:C/H-only','C:O-attached','C:N-attached','H:C-bound','H:O-bound','H:N-bound','O:no-metal','O:metal-adjacent','N:no-metal','N:metal-adjacent']
    fig,axs=plt.subplots(1,3,figsize=(7.4,5.0),sharey=True)
    for j,(q,title,cmap) in enumerate([('calibrated_signed','Signed departure','RdBu_r'),('calibrated_absolute','Absolute departure','YlGnBu'),('residual_sd','Within-environment spread','YlGnBu')]):
        ax=axs[j];mat=atlas[atlas.environment.isin(envs)].pivot(index='environment',columns='pair',values='mean_structure_'+q).reindex(index=envs,columns=PAIRS)
        lim=np.abs(mat.to_numpy()).max()
        im=ax.pcolormesh(np.arange(4)-.5,np.arange(len(envs)+1)-.5,mat.to_numpy(),cmap=cmap,vmin=-lim if j==0 else 0,vmax=lim,shading='flat')
        ax.set_ylim(len(envs)-.5,-.5);ax.set_xlim(-.5,2.5)
        for y in range(len(envs)):
            for x in range(3):
                v=mat.iloc[y,x];ax.text(x,y,f'{v:.3f}',ha='center',va='center',fontsize=7,color='white' if (abs(v)>.6*lim if j==0 else v>.5*lim) else '#20262C')
        ax.set(title=title,xticks=range(3),xticklabels=['DDEC6–CM5','DDEC6–Bader','CM5–Bader'])
        ax.tick_params(axis='x',rotation=50)
        if j==0:ax.set_yticks(range(len(envs)),envs)
        cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.05,pad=.16,label='e')
        cb.solids.set_rasterized(False);panel(ax,chr(97+j))
    save(fig,'Figure_S04_atlas')
    d=read('nonmetal_contrasts.csv')
    fig,axs=plt.subplots(2,3,figsize=(7.4,6.4),sharey=True)
    for i,q in enumerate(['calibrated_signed','calibrated_absolute']):
        for j,pair in enumerate(PAIRS):
            ax=axs[i,j]
            for k,(scale,color,marker) in enumerate([(1.2,BLUE,'v'),(1.25,TEAL,'o'),(1.3,ORANGE,'s')]):
                for y,left in enumerate(LEFT):
                    r=d[d.scale.eq(scale)&d.pair.eq(pair)&d.left.eq(left)&d.quantity.eq(q)].iloc[0]
                    point(ax,r.effect,r.effect_ci_low,r.effect_ci_high,y+(k-1)*.17,color,marker)
            ax.axvline(0,color=GRAY,lw=.7);ax.set(title=NAMES[j],xlabel=('Signed' if i==0 else 'Absolute')+' contrast (e)')
            if j==0:
                ns=d[d.scale.eq(1.25)&d.pair.eq(pair)&d.quantity.eq(q)].set_index('left').n_common_structures if 'n_common_structures' in d else None
                ax.set_yticks(range(6), [label+(f' (n={int(ns[left])})' if ns is not None else '') for label,left in zip(LABELS,LEFT)])
            panel(ax,chr(97+i*3+j))
    axs[0,0].invert_yaxis()
    fig.legend(loc='outside lower center',ncol=3,handles=[Line2D([],[],color=c,marker=m,linestyle='none',label=f'{s:.2f}'+(' stress' if s==1.3 else '')) for s,c,m in [(1.2,BLUE,'v'),(1.25,TEAL,'o'),(1.3,ORANGE,'s')]],frameon=False)
    save(fig,'Figure_S05_cutoff_contrasts')

def main():
    fig1();fig2();contrast_grid(read('contrasts.csv'),'Figure_03_chemical_contrasts');fig4();fig5();supplement()
    (OUT/'figure_build.json').write_text(json.dumps({'inputs':sorted(USED),'figures':sorted(p.stem for p in OUT.glob('*.pdf')),
        'scope':'Frozen saved estimates only; no refit or new bootstrap; table intervals preserved.'},indent=2)+'\n')
    exports=[dict(file=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir()) if p.suffix in ['.pdf','.svg','.png']]
    (ROOT/'figures/FIGURE_OUTPUTS.json').write_text(json.dumps(exports,indent=2)+'\n',encoding='utf-8')
    (ROOT/'figures/figure_environment.json').write_text(json.dumps(dict(python=sys.version.split()[0],packages={n:importlib.metadata.version(n) for n in ['numpy','pandas','matplotlib']},note='Actual figure-build environment; distinct from original fitted-model runtime.'),indent=2)+'\n',encoding='utf-8')
    print(f'Built {len(list(OUT.glob("*.pdf")))} PDF/SVG/PNG figure triplets in {OUT}')

if __name__=='__main__':main()
