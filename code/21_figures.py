"""Generate review figures from saved numerical tables; save exact plotting inputs."""
import os
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';IN=OUT/'inputs'
PAIRS=['ddec6_cm5','ddec6_bader','cm5_bader'];TITLES=['DDEC6 − CM5','DDEC6 − Bader','CM5 − Bader']
COLORS=['#718096','#e69f00','#009e73','#0072b2','#cc79a7']

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight',facecolor='white');fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)

def main():
    if (ROOT/'results/closeout_v1').exists():
        from importlib import import_module
        return import_module('27_closeout_review').make_figures()
    OUT.mkdir(exist_ok=True);IN.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.dpi':110})
    ed=pd.read_csv(ROOT/'results/charge/element_distributions.csv');ss=pd.read_csv(ROOT/'results/charge/strata_summary.csv')
    elements=['H','C','N','O','Mg','Al','Ca','Mn','Fe','Co','Ni','Cu','Zn','Zr','Cd','Pb']
    ed=ed[ed.element.isin(elements)];ed.to_csv(IN/'element_differences.csv',index=False);ss.to_csv(IN/'sign_opposition.csv',index=False)
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7),gridspec_kw={'width_ratios':[1.3,1]})
    mat=ed.pivot(index='element',columns='pair',values='mean_signed').reindex(elements)[PAIRS]
    lim=max(abs(mat.min().min()),abs(mat.max().max()));im=ax.imshow(mat,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
    ns=ed.drop_duplicates('element').set_index('element')
    ax.set_yticks(range(len(elements)),[f'{e}  ({ns.loc[e,"n_structures"]:,} structures)' for e in elements]);ax.set_xticks(range(3),TITLES,rotation=20,ha='right')
    for i in range(len(elements)):
        for j in range(3):
            v=mat.iloc[i,j];ax.text(j,i,f'{v:+.2f}',ha='center',va='center',color='white' if abs(v)>.65*lim else 'black',fontsize=9)
    ax.set_title('Mean signed scheme difference (e)');fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
    x=np.arange(4)
    for k,label in enumerate(['all','nonmetal','metal']):
        row=ss[ss.stratum==label].iloc[0];values=[row[f'sign_opposition_{e:.2f}']*100 for e in (0.,.01,.05,.10)]
        bx.bar(x+(k-1)*.24,values,width=.23,label=f'{label} (n={row.n_atoms:,})',color=COLORS[k])
    bx.set_xticks(x,['0','0.01','0.05','0.10']);bx.set_xlabel('Required magnitude on both sides of zero (e)');bx.set_ylabel('Atoms with opposite-sign charges (%)');bx.set_title('Near-zero sign flips change the count');bx.legend(fontsize=8);bx.set_ylim(0,17)
    fig.tight_layout();save(fig,'01_scheme_differences')

    cr=pd.read_csv(ROOT/'results/charge/calibration_residual_deciles.csv');cr=cr[cr.element.isin(['C','O'])];cr.to_csv(IN/'calibration_curvature.csv',index=False)
    fig,axs=plt.subplots(2,3,figsize=(12,6),sharey='row')
    for i,el in enumerate(['C','O']):
        for j,pair in enumerate(PAIRS):
            g=cr[(cr.element==el)&(cr.pair==pair)].sort_values('qA_mean');ax=axs[i,j]
            ax.plot(g.qA_mean,g.affine_signed_residual,'o-',color=COLORS[1],label='Affine');ax.plot(g.qA_mean,g.nonlinear_signed_residual,'s-',color=COLORS[2],label='12-bin piecewise')
            ax.axhline(0,color='k',lw=.6);ax.set_xlabel(pair.split('_')[0].upper()+' charge (e)');ax.set_title(f'{el}: {TITLES[j]}')
            if j==0:ax.set_ylabel('Mean held-out signed residual (e)')
            if i==0 and j==2:ax.legend(fontsize=9)
    fig.tight_layout();save(fig,'02_calibration_residuals')

    cc=pd.read_csv(ROOT/'results/science/paired_environment_contrasts.csv');cc=cc[(cc.quantity=='signed')&cc.control.isin(['affine','nonlinear'])&cc.supported_denominator];cc.to_csv(IN/'chemical_contrasts.csv',index=False)
    labels=[('C:O-attached','C:C/H-only'),('C:N-attached','C:C/H-only'),('H:O-bound','H:C-bound'),('H:N-bound','H:C-bound'),('O:metal-adjacent','O:no-metal'),('N:metal-adjacent','N:no-metal')]
    fig,axs=plt.subplots(1,3,figsize=(14,5),sharey=True)
    for j,pair in enumerate(PAIRS):
        ax=axs[j]
        for k,control in enumerate(['affine','nonlinear']):
            for i,(left,right) in enumerate(labels):
                gg=cc[(cc.pair==pair)&(cc.control==control)&(cc.left==left)&(cc.right==right)]
                if gg.empty:continue
                r=gg.iloc[0];ax.errorbar(r.effect_e,i+(k-.5)*.18,xerr=[[r.effect_e-r.ci_low],[r.ci_high-r.effect_e]],fmt='o',ms=4,capsize=2,color=COLORS[k+1],label=control if i==0 else None)
        ax.axvline(0,color='k',lw=.7);ax.set_title(TITLES[j]);ax.set_xlabel('Within-structure signed residual contrast (e)')
        if j==0:
            display=[]
            for left,right in labels:
                rr=cc[(cc.pair==pair)&(cc.control=='nonlinear')&(cc.left==left)&(cc.right==right)]
                n=int(rr.iloc[0].n_paired_structures) if len(rr) else 0;display.append(left+' vs '+right.split(':')[1]+f'\n(n={n:,} paired structures)')
            ax.set_yticks(range(len(labels)),display,fontsize=9);ax.invert_yaxis()
        if j==2:ax.legend()
    fig.tight_layout();save(fig,'03_environment_contrasts')

    mm=pd.read_csv(ROOT/'results/science/grouped_metrics.csv');mods=['offset','affine','nonlinear','affine_ridge_first','affine_tree_second','nonlinear_tree_first','nonlinear_tree_second']
    mm=mm[mm.model.isin(mods)&mm.stratum.isin(['all','metal'])];mm.to_csv(IN/'predictive_metrics.csv',index=False)
    domain=pd.read_csv(ROOT/'results/metal_control/grouped_metrics.csv');domain.to_csv(IN/'metal_domain_metrics.csv',index=False)
    fig,axs=plt.subplots(2,3,figsize=(13,8))
    short=['Element','Affine','Piecewise','Affine +\nridge first','Affine +\ntree two-hop','Piecewise +\ntree first','Piecewise +\ntree two-hop']
    for i,stratum in enumerate(['all','metal']):
        for j,pair in enumerate(PAIRS):
            ax=axs[i,j];g=mm[(mm.pair==pair)&(mm.stratum==stratum)].set_index('model').reindex(mods)
            if stratum=='all':
                ax.errorbar(range(len(mods)),g.mae_structure,yerr=[g.mae_structure-g.mae_structure_ci_low,g.mae_structure_ci_high-g.mae_structure],fmt='o-',color=COLORS[3],ms=4,capsize=3)
                ax.set_xticks(range(len(mods)),short,rotation=65,ha='right',fontsize=8)
            else:
                dg=domain[(domain.pair==pair)&domain.stratum.eq('all_metals')].set_index('model')
                values=[g.loc['nonlinear'],dg.loc['element_residual'],g.loc['nonlinear_tree_first'],g.loc['nonlinear_tree_second'],dg.loc['metal_tree_first'],dg.loc['metal_tree_second']]
                for k,r in enumerate(values):
                    low=r['mae_structure_ci_low'] if k in (0,2,3) else r['mae_ci_low'];high=r['mae_structure_ci_high'] if k in (0,2,3) else r['mae_ci_high']
                    ax.errorbar(k,r.mae_structure,yerr=[[r.mae_structure-low],[high-r.mae_structure]],fmt='o',color=COLORS[3] if k in (0,2,3) else COLORS[2],ms=5,capsize=3)
                ax.set_xticks(range(6),['Piecewise','Metal element\nresidual offset','Global tree\nfirst','Global tree\ntwo-hop','Metal tree\nfirst','Metal tree\ntwo-hop'],rotation=65,ha='right',fontsize=8)
            ax.set_title(TITLES[j]+' | '+stratum);ax.set_ylim(bottom=0)
            if j==0:ax.set_ylabel('Structure-balanced held-out MAE (e)')
    fig.tight_layout();save(fig,'04_grouped_prediction')

    sens=pd.read_csv(ROOT/'results/sensitivity/heldout_graph_metrics.csv');sens=sens[(sens.model=='nonlinear_tree_second')&sens.stratum.eq('all')]
    ds=pd.read_csv(ROOT/'results/metal_control/graph_sensitivity.csv');ds=ds[ds.model.eq('metal_tree_second')].rename(columns={'ci_low':'delta_ci_low','ci_high':'delta_ci_high'});ds['stratum']='metal'
    sens=pd.concat([sens,ds],ignore_index=True);sens.to_csv(IN/'graph_sensitivity_metrics.csv',index=False)
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1,1.3]})
    records={r['qmof_id']:r for r in json.loads((ROOT/'inputs/expanded_benchmark_records.json').read_text())};r=records['qmof-8a95c27'];from chamm_graph import structure
    s=structure(r);edges=json.loads((ROOT/'results/connectivity/qmof-8a95c27.json').read_text())['methods']['hybrid_metal_center_union'][2]
    v=np.array([(s[e[0]].frac_coords+np.array(e[1:4])-s[2].frac_coords)@s.lattice.matrix for e in edges]);_,_,vh=np.linalg.svd(v);xy=v@vh[:2].T
    for p,e in zip(xy,edges):
        col=COLORS[1] if e[4]>2.3 else COLORS[3];ax.plot([0,p[0]],[0,p[1]],color=col);ax.scatter(*p,s=70,c=col,zorder=3)
        offset=(6,-12) if p[0]<-1.5 else ((-5,5) if p[0]<0 else (5,5))
        ax.annotate(f'O {e[4]:.3f} Å',p,xytext=offset,ha='left' if p[0]<-1.5 or p[0]>=0 else 'right',textcoords='offset points',fontsize=8)
    ax.scatter(0,0,s=140,c='#8c564b',zorder=4);ax.annotate('Cu',(0,0),xytext=(5,-15),textcoords='offset points');ax.set_aspect('equal');ax.set_xlabel('Local principal projection 1 (Å)');ax.set_ylabel('Projection 2 (Å)');ax.set_title('Cu2: plausible longer opposite pair\n2.3985 Å, angle 179.999°')
    yy=np.arange(6)
    for k,variant in enumerate(['shift','cutoff']):
        for j,pair in enumerate(PAIRS):
            for i,stratum in enumerate(['all','metal']):
                row=sens[(sens.variant==variant)&(sens.pair==pair)&(sens.stratum==stratum)].iloc[0]
                y=j*2+i+(k-.5)*.18
                bx.errorbar(row.alternative_minus_accepted,y,xerr=[[row.alternative_minus_accepted-row.delta_ci_low],[row.delta_ci_high-row.alternative_minus_accepted]],fmt='o',color=COLORS[k+1],capsize=3,label=('Taper +0.10 Å' if variant=='shift' else 'Original cutoff stress') if j==0 and i==0 else None)
    bx.axvline(0,color='k',lw=.7);bx.set_yticks(yy,[f'{TITLES[j]} | {st}' for j in range(3) for st in ['all (global)','metals (domain)']],fontsize=9);bx.invert_yaxis();bx.set_xlabel('Alternative − accepted held-out MAE (e)');bx.set_title('96 structures: fixed models, changed graph');bx.legend(fontsize=8)
    fig.tight_layout();save(fig,'05_connectivity_sensitivity')

if __name__=='__main__':main()
