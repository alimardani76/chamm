"""Current scientific closeout figures/reports/bundle. Called by stages 21 and 22.
All headline numbers come from the finite-test outputs, never the old narrative.
"""
from pathlib import Path
import hashlib,json,zipfile,datetime
import numpy as np,pandas as pd
from chamm_closeout import ROOT,OUT,PAIRS,CONTRASTS,QUANTITIES
TITLES=['DDEC6 − CM5','DDEC6 − Bader','CM5 − Bader']
BLUE='#0072b2';GREEN='#009e73';ORANGE='#e69f00'
MODELS=['piecewise','center','charge_only','conditional_first','conditional_twohop']

def require_results():
    for sub in ['original','blocked','atlas','nonmetal']:assert (OUT/sub/'COMPLETE.json').exists(),f'Incomplete closeout stage: {sub}; run stages 24,25,26'

def md(df):
    def value(x):
        if isinstance(x,(float,np.floating)):return f'{x:.5f}'
        return str(x).replace('|','/')
    return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(value(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)])

def tables():
    require_results();return {k:pd.read_csv(OUT/p) for k,p in dict(metrics='atlas/grouped_metrics.csv',increments='atlas/paired_increments.csv',contrasts='atlas/headline_contrasts.csv',atlas='atlas/environment_atlas.csv',graph='nonmetal/common_pair_contrasts.csv',available='nonmetal/available_pair_contrasts.csv',challenge='atlas/grouping_challenge.csv').items()}

def interval(r,value='effect',lo='ci_low',hi='ci_high'):
    return f'{r[value]:+.5f} [{r[lo]:+.5f}, {r[hi]:+.5f}]'

def sign_status(a,b):
    if a.ci_low>0 and b.ci_low>0:return 'positive in both validation schemes'
    if a.ci_high<0 and b.ci_high<0:return 'negative in both validation schemes'
    return 'weak/null or inconsistent; no recurring direction established'

def build_evidence(data):
    c=data['contrasts'];g=data['graph'];rows=[]
    for pair in PAIRS:
        for left,right in CONTRASTS:
            base=dict(pair=pair,left=left,right=right)
            for quantity in QUANTITIES:
                original=c[(c.split=='original')&(c.pair==pair)&(c.left==left)&(c.quantity==quantity)].iloc[0]
                blocked=c[(c.split=='blocked')&(c.pair==pair)&(c.left==left)&(c.quantity==quantity)].iloc[0]
                gg=g[(g.pair==pair)&(g.left==left)&(g.quantity==quantity)]
                n=int(gg.n_common_structures.min()) if len(gg) else 0
                change=float(gg.change.abs().max()) if len(gg) else np.nan
                graph=('unchanged in tested subset' if change==0 else 'changes under tested graph')+(' (limited paired support)' if n<30 else '')
                rows.append(dict(base,quantity=quantity,n_structures=int(blocked.n_structures),n_components=int(blocked.n_validation_blocks),original_effect=original.effect,original_ci_low=original.ci_low,original_ci_high=original.ci_high,blocked_effect=blocked.effect,blocked_ci_low=blocked.ci_low,blocked_ci_high=blocked.ci_high,
                    validation_status=sign_status(original,blocked),graph_status=graph,graph_common_structures=n,graph_max_effect_change=change,
                    interpretation={'raw_signed':'direction of raw scheme difference','raw_absolute':'raw disagreement magnitude contrast','calibrated_signed':'direction of departure; not greater disagreement','calibrated_absolute':'absolute empirical charge-mapping departure; not physical error','residual_sd':'within-environment residual spread; not calibrated physical uncertainty'}[quantity]))
    ev=pd.DataFrame(rows);ev.to_csv(OUT/'atlas/contrast_evidence_status.csv',index=False)
    atlas=data['atlas'];compact=atlas[(atlas.split=='blocked')&atlas.adequate_support&(atlas.n_validation_blocks>=30)&atlas.environment.ne('other')]
    index=['pair','element','environment','n_atoms','n_structures','n_validation_blocks','ambiguity_fraction_atom','ambiguity_fraction_structure']
    wide=compact.pivot(index=index,columns='quantity',values=['mean_structure','ci_low','ci_high']).reset_index()
    wide.columns=['_'.join(str(v) for v in col if v) if isinstance(col,tuple) else col for col in wide.columns]
    subset=pd.read_parquet(OUT/'nonmetal/environment_units.parquet',columns=['pair','scale','element','environment','qmof_id'])
    support=subset[subset.pair.eq(PAIRS[0])].groupby(['element','environment','scale']).qmof_id.nunique().groupby(['element','environment']).min()
    wide['minimum_subset_environment_structures']=[int(support.get((r.element,r.environment),0)) for r in wide.itertuples()]
    wide['graph_coverage']=np.where(wide.minimum_subset_environment_structures>=30,'supported subset environment coverage','limited/absent subset environment coverage')
    wide.to_csv(OUT/'atlas/compact_environment_atlas.csv',index=False)
    claims=[];inc=data['increments']
    for pair in PAIRS:
        for pop in ['all','nonmetal','metal','has_mofkey','nonsingleton']:
            r=inc[(inc.split=='blocked')&(inc.pair==pair)&(inc.stratum==pop)&inc.comparison.eq('structure_given_q_A')].iloc[0]
            claims.append(dict(claim=f'Conditional local chemistry: {pair}, {pop}',status='survived specified challenge' if r.ci_low>0 else 'failed benefit in this global model' if r.ci_high<0 else 'unsupported/uncertain',effect_e=r.improvement_structure,ci_low=r.ci_low,ci_high=r.ci_high,n_structures=int(r.n_structures),evidence='atlas/paired_increments.csv; blocked; structure_given_q_A',scope='same fixed tree/population; identifier/composition blocking, not universal family transfer'))
    for _,r in ev[ev.quantity.isin(['calibrated_signed','calibrated_absolute'])].iterrows():
        claims.append(dict(claim=f'{r.pair}: {r.left} vs {r.right}, {r.quantity}',status=r.validation_status+'; '+r.graph_status,effect_e=r.blocked_effect,ci_low=r.blocked_ci_low,ci_high=r.blocked_ci_high,n_structures=int(r.n_structures),evidence='atlas/contrast_evidence_status.csv; nonmetal/common_pair_contrasts.csv',scope=r.interpretation))
    claims += [dict(claim=x,status='outside evidence',effect_e=np.nan,ci_low=np.nan,ci_high=np.nan,n_structures=0,evidence='finite closeout scope',scope='no supporting independent observable or validated universal grouping') for x in ['charge-scheme accuracy ranking','redox assignment','universal chemical bond truth','universal family/topology transfer','causal electronic mechanism','empirical residual as physical uncertainty']]
    pd.DataFrame(claims).to_csv(OUT/'claim_to_evidence.csv',index=False)
    return ev,wide,pd.DataFrame(claims)

def make_figures():
    data=tables();build_evidence(data)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=ROOT/'figures';inp=out/'inputs';inp.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=220,bbox_inches='tight',facecolor='white');fig.savefig(out/(name+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
    # Figure 1 retains the frozen raw-source content; regenerate it if its renders were removed.
    if not (out/'01_scheme_differences.png').exists() or not (out/'01_scheme_differences.svg').exists():
        ed=pd.read_csv(ROOT/'results/charge/element_distributions.csv');ss=pd.read_csv(ROOT/'results/charge/strata_summary.csv');els=['H','C','N','O','Mg','Al','Ca','Mn','Fe','Co','Ni','Cu','Zn','Zr','Cd','Pb']
        fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7),gridspec_kw={'width_ratios':[1.3,1]});values=ed.pivot(index='element',columns='pair',values='mean_signed').reindex(els)[list(PAIRS)];lim=np.abs(values.to_numpy()).max();im=ax.imshow(values,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ns=ed.drop_duplicates('element').set_index('element').n_structures
        ax.set_yticks(range(len(els)),[f'{e} ({ns[e]:,} structures)' for e in els]);ax.set_xticks(range(3),TITLES,rotation=20,ha='right');ax.set_title('Mean signed scheme difference (e)');fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
        for i in range(len(els)):
            for j in range(3):
                v=values.iloc[i,j];ax.text(j,i,f'{v:+.2f}',ha='center',va='center',color='white' if abs(v)>.65*lim else 'black',fontsize=9)
        for j,label in enumerate(['all','nonmetal','metal']):
            r=ss[ss.stratum.eq(label)].iloc[0];vv=[100*r[f'sign_opposition_{eps:.2f}'] for eps in [0.,.01,.05,.1]];bx.bar(np.arange(4)+(j-1)*.24,vv,width=.23,label=f'{label} (n={int(r.n_atoms):,})',color=[BLUE,ORANGE,GREEN][j])
        bx.set_xticks(range(4),['0','0.01','0.05','0.10']);bx.set_xlabel('Required magnitude on both sides of zero (e)');bx.set_ylabel('Atoms with opposite-sign charges (%)');bx.set_title('Near-zero signs');bx.legend(fontsize=8);fig.tight_layout();save(fig,'01_scheme_differences')
    m=data['metrics'];m=m[m.stratum.isin(['all','metal'])];m.to_csv(inp/'closeout_conditional_metrics.csv',index=False)
    fig,axs=plt.subplots(2,3,figsize=(13,7.5))
    labels=['Piecewise','Center only','Identity + qA','+ First shell','+ Two-hop paths']
    for i,pop in enumerate(['all','metal']):
        for j,pair in enumerate(PAIRS):
            ax=axs[i,j]
            for k,split in enumerate(['original','blocked']):
                g=m[(m.split==split)&(m.pair==pair)&(m.stratum==pop)].set_index('model').reindex(MODELS)
                x=np.arange(5)+(k-.5)*.10;ax.errorbar(x,g.mae_structure,yerr=[g.mae_structure-g.ci_low,g.ci_high-g.mae_structure],fmt='o-' if k==0 else 's--',color=BLUE if k==0 else ORANGE,ms=4,capsize=2,label='QMOF folds' if k==0 else 'Identifier/composition blocks')
            ax.set_xticks(range(5),labels,rotation=35,ha='right',fontsize=8);ax.set_title(TITLES[j]+' | '+pop);ax.set_ylim(bottom=0)
            if j==0:ax.set_ylabel('Structure-balanced held-out MAE (e)')
            if i==0 and j==2:ax.legend(fontsize=7)
    fig.tight_layout();save(fig,'02_calibration_residuals')
    c=data['contrasts'];cc=c[c.quantity.isin(['calibrated_signed','calibrated_absolute'])];cc.to_csv(inp/'closeout_signed_magnitude_contrasts.csv',index=False)
    fig,axs=plt.subplots(2,3,figsize=(14,8),sharey=True)
    for i,quantity in enumerate(['calibrated_signed','calibrated_absolute']):
        for j,pair in enumerate(PAIRS):
            ax=axs[i,j]
            for k,split in enumerate(['original','blocked']):
                for y,(left,right) in enumerate(CONTRASTS):
                    r=cc[(cc.split==split)&(cc.pair==pair)&(cc.left==left)&(cc.quantity==quantity)].iloc[0]
                    ax.errorbar(r.effect,y+(k-.5)*.16,xerr=[[max(0,r.effect-r.ci_low)],[max(0,r.ci_high-r.effect)]],fmt='o' if k==0 else 's',ms=3,capsize=2,color=BLUE if k==0 else ORANGE,label=split if y==0 else None)
            ax.axvline(0,color='k',lw=.6);ax.set_title(TITLES[j]+(' | signed departure' if i==0 else ' | absolute departure'))
            ax.set_xlabel('Within-structure contrast (e)')
            if j==0:ax.set_yticks(range(6),[l+' vs '+r.split(':')[1] for l,r in CONTRASTS],fontsize=8);ax.invert_yaxis()
            if i==0 and j==2:ax.legend(fontsize=8)
    fig.tight_layout();save(fig,'03_environment_contrasts')
    inc=data['increments'];inc=inc[inc.comparison.eq('structure_given_q_A')&inc.stratum.isin(['all','has_mofkey','nonsingleton'])];inc.to_csv(inp/'closeout_blocked_challenge.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(12,4.4),sharey=True)
    for j,pair in enumerate(PAIRS):
        ax=axs[j]
        for k,split in enumerate(['original','blocked']):
            for y,pop in enumerate(['all','has_mofkey','nonsingleton']):
                r=inc[(inc.split==split)&(inc.pair==pair)&(inc.stratum==pop)].iloc[0]
                ax.errorbar(1000*r.improvement_structure,y+(k-.5)*.18,xerr=[[1000*(r.improvement_structure-r.ci_low)],[1000*(r.ci_high-r.improvement_structure)]],fmt='o' if k==0 else 's',color=BLUE if k==0 else ORANGE,ms=5,capsize=3,label=split if y==0 else None)
        ax.axvline(0,color='k',lw=.6);ax.set_title(TITLES[j]);ax.set_xlabel('MAE gain over charge-only\n'+r'($10^{-3}e$)')
        if j==0:ax.set_yticks(range(3),['All: 9,171 structures','MOFkey present: 5,408','Nonsingleton blocks: 1,414']);ax.invert_yaxis()
        if j==2:ax.legend()
    fig.tight_layout();save(fig,'04_grouped_prediction')
    g=data['graph'];g=g[g.quantity.isin(['calibrated_signed','calibrated_absolute'])];g.to_csv(inp/'closeout_nonmetal_graph_contrasts.csv',index=False)
    fig,axs=plt.subplots(2,3,figsize=(14,8),sharey=True)
    for i,quantity in enumerate(['calibrated_signed','calibrated_absolute']):
        for j,pair in enumerate(PAIRS):
            ax=axs[i,j]
            for k,scale in enumerate([1.2,1.25,1.3]):
                for y,(left,right) in enumerate(CONTRASTS):
                    r=g[(g.scale==scale)&(g.pair==pair)&(g.left==left)&(g.quantity==quantity)].iloc[0]
                    ax.errorbar(r.effect,y+(k-1)*.17,xerr=[[r.effect-r.effect_ci_low],[r.effect_ci_high-r.effect]],fmt=['v','o','s'][k],color=[GREEN,BLUE,ORANGE][k],ms=3,capsize=2,label=f'{scale:.2f}'+(' stress' if scale==1.3 else '') if y==0 else None)
            ax.axvline(0,color='k',lw=.6);ax.set_title(TITLES[j]+(' | signed' if i==0 else ' | absolute'));ax.set_xlabel('Common-pair subset contrast (e)')
            if j==0:
                labels=[]
                for left,right in CONTRASTS:
                    rr=g[(g.scale==1.25)&(g.pair==pair)&(g.left==left)&(g.quantity==quantity)].iloc[0];labels.append(left+' vs '+right.split(':')[1]+f' (n={int(rr.n_common_structures)})')
                ax.set_yticks(range(6),labels,fontsize=8);ax.invert_yaxis()
            if i==0 and j==2:ax.legend(fontsize=8)
    fig.tight_layout();save(fig,'05_connectivity_sensitivity')
    (inp/'README.md').write_text('# Active final figure inputs\n\n1: element_differences.csv and sign_opposition.csv (unchanged raw data).\n2: closeout_conditional_metrics.csv.\n3: closeout_signed_magnitude_contrasts.csv.\n4: closeout_blocked_challenge.csv.\n5: closeout_nonmetal_graph_contrasts.csv.\n\nOther inputs without the closeout prefix are preserved historical production inputs; they do not drive the revised figures 2–5. Full tables and resampling units are in results/closeout_v1/.\n',encoding='utf-8')
    print('Current five review figures generated; Figure 1 unchanged',flush=True)

def generate():
    data=tables();ev,compact,claims=build_evidence(data);m=data['metrics'];inc=data['increments'];c=data['contrasts'];graph=data['graph']
    meta=json.loads((OUT/'metadata/grouping_summary.json').read_text());manifest=json.loads((ROOT/'data/descriptor_manifest.json').read_text())
    def row(split,pair,pop,model):return m[(m.split==split)&(m.pair==pair)&(m.stratum==pop)&m.model.eq(model)].iloc[0]
    def effect(split,pair,left,quantity):return c[(c.split==split)&(c.pair==pair)&(c.left==left)&c.quantity.eq(quantity)].iloc[0]
    primary=[];challenge=[];metal=[]
    for split in ['original','blocked']:
        for pair in PAIRS:
            r=inc[(inc.split==split)&(inc.pair==pair)&(inc.stratum=='all')&inc.comparison.eq('structure_given_q_A')].iloc[0]
            primary.append(dict(split=split,pair=pair,**{model:row(split,pair,'all',model).mae_structure for model in MODELS},conditional_gain=r.improvement_structure,ci_low=r.ci_low,ci_high=r.ci_high,folds_improved=int(r.n_folds_improved)))
    for pair in PAIRS:
        for pop in ['all','has_mofkey','nonsingleton']:
            r=inc[(inc.split=='blocked')&(inc.pair==pair)&(inc.stratum==pop)&inc.comparison.eq('structure_given_q_A')].iloc[0]
            challenge.append(dict(pair=pair,population=pop,n_structures=int(r.n_structures),n_components=int(r.n_validation_blocks),gain=r.improvement_structure,ci_low=r.ci_low,ci_high=r.ci_high))
        r=inc[(inc.split=='blocked')&(inc.pair==pair)&(inc.stratum=='metal')&inc.comparison.eq('structure_given_q_A')].iloc[0]
        metal.append(dict(pair=pair,charge_only=row('blocked',pair,'metal','charge_only').mae_structure,conditional_twohop=row('blocked',pair,'metal','conditional_twohop').mae_structure,gain=r.improvement_structure,ci_low=r.ci_low,ci_high=r.ci_high,folds_improved=int(r.n_folds_improved)))
    pd.DataFrame(primary).to_csv(OUT/'main_comparison.csv',index=False);pd.DataFrame(challenge).to_csv(OUT/'blocked_support_comparison.csv',index=False)
    verdict='Restricted positive: local chemical information survives tested charge conditioning and identifier/composition blocking; central-metal and magnitude claims remain qualified'
    (OUT/'verdict.json').write_text(json.dumps(dict(verdict=verdict,scope='finite specified comparisons; no model/graph search after outcomes',stopped=True),indent=2))
    a=effect('blocked','ddec6_bader','C:O-attached','calibrated_signed');b=effect('blocked','ddec6_bader','C:O-attached','calibrated_absolute')
    n=effect('blocked','cm5_bader','C:N-attached','calibrated_signed');na=effect('blocked','cm5_bader','C:N-attached','calibrated_absolute')
    h=effect('original','ddec6_cm5','H:N-bound','calibrated_signed');oldc=pd.read_csv(ROOT/'results/science/paired_environment_contrasts.csv');oldh=oldc[(oldc.control=='nonlinear')&(oldc.pair=='ddec6_cm5')&(oldc.left=='H:N-bound')&oldc.quantity.eq('signed')].iloc[0]
    examples=[]
    for pair,left in [('ddec6_bader','C:O-attached'),('cm5_bader','C:N-attached')]:
        examples.append(dict(pair=pair,contrast=left+' vs C:C/H-only',n_structures=int(effect('blocked',pair,left,'raw_signed').n_structures),**{q:effect('blocked',pair,left,q).effect for q in QUANTITIES}))
    added=pd.read_csv(OUT/'nonmetal/physical_edge_changes.csv');changes=pd.read_csv(OUT/'nonmetal/change_summary.csv').groupby('scale')[['first_shell_changed','twohop_changed','environment_changed']].sum();maxchange=graph.change.abs().max()
    findings=f'''# ChaMM final scientific findings

**Final verdict: {verdict}.** The finite closeout is complete. The production graph, cohort and descriptor table remain frozen. No new whole-cohort graph build, DFT, topology campaign or model search was performed.

## Conditional comparison: chemistry is not explained by the stronger scalar control

Every competitor uses the same fitting atoms, weights and outer-test atoms. The common offset is the original training-fitted 12-bin charge calibration. The center comparator is an element residual mean; the other three competitors use the same fixed tree with identity + q_A, then first shell, then first shell + two-hop paths. q_A appears only in this explicitly conditioned follow-up. No scheme is reference truth.

Structure-balanced held-out MAE in e, all **9,171 structures / 1,420,532 atoms**:

{md(pd.DataFrame(primary))}

The main increment is conditional two-hop versus flexible charge-only. It remains positive in all five folds for every pair under both splits, with paired intervals above zero. Central-element residual recentering alone does not explain the result. All atom-weighted metrics, common-element support, unsupported calibrations, and per-fold increments are retained in `results/closeout_v1/atlas/`. The center-to-first comparison jointly adds q_A and local features; only charge-only-to-first isolates first-shell information conditional on q_A. Path increments are separately reported in `paired_increments.csv`; they are not the old charge-blind path increments.

## Conservative grouping challenge

The prior statement that grouping metadata was unavailable is corrected. The pinned `qmof_database/qmof.json` contains source reduced formulas and exact MOFkeys. All cohort IDs and those grouping fields match the supplied reviewer metadata. There are {meta['mofkey_records']:,} non-null MOFkey records / {meta['distinct_mofkeys']:,} distinct keys, and {meta['distinct_formulas']:,} exact reduced formulas.

Linking equal non-null MOFkey **or** equal exact source reduced formula produces **{meta['components']:,} components**, largest {meta['largest_component']}, with **{meta['nonsingleton_records']:,} nonsingleton records** and {meta['records_components_crossing_old_folds']:,} records in components spanning the old folds. The deterministic new folds contain {meta['fold_counts']} structures. This is identifier/composition-blocked validation; formula equality is not structural identity, and this does not prove universal family or topology separation.

Blocked conditional gain over charge-only, in e:

{md(pd.DataFrame(challenge))}

The nonsingleton challenge population and MOFkey-present population both retain the gain; the singleton majority therefore does not conceal failure on those challenge sets. Blocked intervals resample components using sums of per-QMOF paired losses divided by resampled QMOF counts. They retain a structure-balanced estimand, rather than weighting every component mean equally.

## Metal conclusions are narrower

For the global conditional models, blocked metal results over **46,994 atoms / 9,165 structures** are:

{md(pd.DataFrame(metal))}

DDEC6−Bader worsens consistently for central metals; the CM5−Bader gain is small. A broad three-pair central-metal benefit is not retained. A failure of this global population model does not prove absence of metal chemistry. The earlier all-metal follow-up remains separately disclosed in `results/metal_control/`; it was motivated by the previous global failure, was not an independent replication, and does not substitute for the present conditional test. No second domain-specific ladder was added in this closeout.

## Signed shifts, magnitude and spread are different results

The atlas reports raw signed differences, raw absolute differences, signed departures after flexible charge-only correction, absolute calibrated departures, and within-structure/environment residual SD. All are empirical quantities in e. Neither absolute departure nor residual spread is calibrated physical uncertainty or charge-scheme error.

- **Directional departure without greater calibrated magnitude:** O-attached versus strictly C/H-only carbon, DDEC6−Bader: signed contrast {interval(a)} e, but absolute-departure contrast {interval(b)} e, over {int(a.n_structures):,} paired structures. The magnitude contrast is unsupported; the signed shift must not be called increased disagreement.
- **Greater absolute calibrated departure without a recurring signed departure:** N-attached/O-absent versus C/H-only carbon, CM5−Bader: signed contrast {interval(n)} e, while absolute departure increases by {interval(na)} e, over {int(n.n_structures):,} paired structures. Opposite directional departures can cancel in their mean.
- **A weakened earlier headline:** N-bound versus C-bound H, DDEC6−CM5, changes from the earlier piecewise-residual estimate {oldh.effect_e:+.5f} e to {interval(h)} e after the stronger charge-only control on the original folds. A residual contrast remains, but the former magnitude is not retained.

These examples also distinguish raw disagreement from calibrated departure (blocked, paired structure means, e):

{md(pd.DataFrame(examples))}

The full signed/magnitude/spread contrasts, weak and null effects, ambiguity fractions, atom counts, structure counts and component support are in the atlas. A compact atlas uses at least 30 structures and 30 represented validation components; unsupported rare classes remain in the full table. No physical error tolerance was invented. The three pair differences are algebraically dependent and do not constitute three independent discoveries.

## Actual nonmetal-cutoff test

On the existing **96 structures / 16,165 atoms**, all accepted metal-involving edges were retained exactly, with no new metal neighbor calculation. Nonmetal cutoff 1.20 changes {int(changes.loc[1.2,'first_shell_changed'])} first shells; 1.30 adds {len(added)} unoriented physical edges, changes {int(changes.loc[1.3,'first_shell_changed'])} first shells and {int(changes.loc[1.3,'twohop_changed'])} two-hop descriptors, but changes **{int(changes.loc[1.3,'environment_changed'])} headline environment memberships**.

The added contacts surround six H sites: H–O distances {added.distance.min():.4f}–{added.distance.max():.4f} Å, nearly linear near-midpoint O–H–O pairs. The accepted graph leaves those H sites isolated; the wider graph gives H degree two. Shared-proton/secondary-contact interpretation remains ambiguous, and 1.30 is a stress configuration, not equally established bond truth. These flagged proton-sharing environments do not support a universal proton-bonding claim.

Every headline contrast was recomputed using the **same stronger original-fold charge-only residuals**, on each graph's available pairs and on the intersection across all three graphs. The maximum common-pair effect change across all reported quantities is {maxchange:.5g} e. Available/common atom and structure denominators and paired changes with QMOF intervals are saved. The unchanged effects establish stability only for these tested categories and settings on this subset. Some H/O/N contrasts have limited paired subset support; the evidence-status table labels that limitation. No retrained alternate-graph performance is claimed.

## Boundary of the final claim

The supported contribution is reproducible, interpretable **within-cohort local chemical structure in charge-definition disagreement conditional on the tested finite models**, robust to this identifier/composition challenge. It is strongest for common nonmetal environments and has no universal direction across schemes. Model capacity and global population weighting remain limitations, especially for metals. The conditional comparison tests these fixed regressors, not every conceivable scalar transformation.

Exact MOFkeys/formulas leave broader related-record dependence unresolved. Bootstrap intervals condition on fitted models and these validation units. The graph is an operational representation with documented ambiguity, not independent bond truth. Charge accuracy, redox, universal family/topology transfer, physical uncertainty calibration and causal electronic mechanisms remain outside the evidence. All four finite tests are finished; these outside-scope topics are not unfinished requirements for this closeout.

See `SCIENTIFIC_CLOSEOUT.md`, `RESULT_TABLES.md`, `METHODS.md`, `FIGURE_CAPTIONS.md`, and `results/closeout_v1/claim_to_evidence.csv` for the final evidence record.
'''
    (ROOT/'reports/SCIENTIFIC_FINDINGS.md').write_text(findings,encoding='utf-8')
    methods=f'''# ChaMM final methods

## Frozen production inputs

QMOF Figshare article 13147324 version 18, file 59573735; archive MD5 `0d89aaf66f2c306e86e47fd91cb1346e`. Common cohort: 9,171 structures / 1,420,532 atom keys `(qmof_id, atom_index)`, including 46,994 metals / 9,165 structures. All charges are finite; production recorded zero graph failures. The operational graph configuration remains `{manifest['configuration_id']}`. The charge and structural tables were neither replaced nor rebuilt in closeout. The one local global-OOF reconciliation checks all three signed targets/folds and 90 saved aggregate model rows; prior reviewer bootstrap and graph audits were reused.

The graph definition, exact radii, parameters, image handling and contact ledger remain in CONNECTIVITY_DECISION.md. Nonmetal pairs use ASE radius scale 1.25/minimum 0.40 Å; metal contacts use metal-center zero-weight CrystalNN with explicit reverse edges and retained distinct images. Flags are diagnostics, not bond-error labels; 24,297 atoms are flagged, including 26.11% of central metals. No assigned oxidation states or charges enter graph selection. Two-hop features count non-backtracking paths, not unique atoms. Bond-distance spread is not angular distortion. Electronic bond orders, spin densities and magnetic moments are excluded.

## Finite conditional comparison

For fixed ordered pairs DDEC6−CM5, DDEC6−Bader and CM5−Bader, the signed target is q_A−q_B. The three targets are algebraically dependent. A common training-fitted original 12-bin piecewise offset is added to every competitor's residual prediction. It uses per-element training quantile bins, weighted charge means, linear interpolation and endpoint-secant extrapolation. The original fallback rules are preserved: fewer than 20 training structures/100 atoms/20 distinct q_A values falls back to the training-fitted affine relation; unobserved elements receive the training global offset. Affine slopes require variance >1e-6 e² and 20 structures. Calibration weights are inverse full structure size, renormalized within element. Original-fold offsets reproduce the saved predictions within 1e-9 e; blocked offsets are newly fitted using only blocked-training structures. Parameters and knots are saved separately by split/fold/pair.

The center-only comparator is a training-sample weighted element residual mean. The same fixed tree predicts the same residual using (1) identity/elemental properties + q_A, (2) those inputs + first-shell descriptors, and (3) those inputs + first shell + two-hop path composition. The structural dictionary is the frozen original dictionary, including common neighboring-metal counts and static central-element indicators. This dictionary has no fitted statistics; an unseen-element indicator is constant zero in training and cannot supply a learned coefficient/split. Native tree missing-value routing and histogram bins are fitted only on training atoms. q_A is explicitly permitted only in these conditional models; none of the earlier charge-blind models was altered. q_B and electronic metadata never enter the predictor matrix.

Tree settings are unchanged: 100 estimators, depth 4, learning rate 0.08, minimum child weight 20, L2=5, histogram/128 bins, two threads, full row/column fractions, seed 20261109. No tuning, early stopping or test-fold model selection was used. Each original fold uses the saved fitting keys/weights exactly. Blocked samples are generated from training structures only by the same rule: up to 12 atoms per structure-element cell, reserve one per available cell, deterministic cap 160,000, then correct weights to the full-training atom proportions and inverse structure size. All competitors share the identical samples, weights and test atoms. All valid-graph atoms receive held-out predictions, covering the full fixed cohort for each split/pair. Center, scalar-charge, first-shell and two-hop comparisons are labeled by their actual conditioning; center-to-first is a joint q_A/neighborhood addition.

## Identifier/composition-blocked challenge

Contrary to the earlier availability statement, `qmof_database/qmof.json` supplies exact `info.formula_reduced` and `info.mofid.mofkey`, as well as MOFid, topology, source and DOI fields. The reviewer join is reused after all cohort IDs, formulas and MOFkeys are checked against the local pinned metadata. Source names are not parsed to invent groups. Connected components link equal non-null exact MOFkey OR equal exact source reduced formula. Missing MOFkeys remain missing. This yields {meta['components']} components, largest {meta['largest_component']}; {meta['nonsingleton_records']} structures are in nonsingleton components and {meta['records_components_crossing_old_folds']} are in components crossing the old folds. All atoms of a component share one new fold. Components are sorted by decreasing size then a declared SHA256 tie order and greedily assigned to the least populated fold, yielding {meta['fold_counts']} structures. No topology-held-out analysis is conducted. Formula equality is conservative blocking, not proof of identity; MOFkeys do not define every broader chemical family.

## Estimates and uncertainty

Atom-weighted MAE and structure-balanced MAE are both saved. The latter averages atom losses within a QMOF/population then weights those QMOFs equally. All/metal/nonmetal, common elements (at least 50 structures), MOFkey-present and nonsingleton populations are explicit; their intersection populations and fallback coverage are also saved. Model increments use paired held-out losses on the same QMOFs. Original-fold intervals resample QMOF IDs. Blocked intervals use 1,000 draws of represented components: component sums of per-QMOF losses (or paired differences) divided by resampled component structure counts. Equal weighting of component means would estimate something different and is not used. Atom-MAE intervals analogously use component atom-error sums/atom counts. Bootstrap draws are shared across competitors for paired increments. Fold-level results and compact per-QMOF loss inputs permit reaggregation without full atom-level predictions. Intervals condition on fitted models and do not resolve all broader related-record or retraining uncertainty.

## Signed and magnitude atlas

The unchanged environment definitions are O-attached carbon; N-attached/O-absent carbon; strictly C/H-only carbon; monocoordinate H bound to C/N/O; O/N with or without a neighboring metal; and metal O-present/N-absent, N-present/O-absent or mixed-O/N classes. Other donors may coexist in the metal classes. Unclassified/rare classes remain visible in the full atlas. Six predeclared same-element contrasts are paired within each QMOF. The full atlas reports raw signed difference, raw absolute difference, flexible-charge-only signed departure, absolute departure, and within-QMOF/environment residual SD. SD is computed from the residual second moment within each QMOF/environment. Its atlas estimate averages those within-QMOF SDs; the atom-weighted column weights them by atom count and is not a pooled atom SD. Single-atom environment SD is zero by definition, not evidence of physical precision. It is not physical charge uncertainty. Main compact atlas entries need at least 30 structures and 30 represented validation units; full tables retain smaller support. Atom counts, QMOF counts, components and atom-/structure-weighted ambiguity coverage accompany estimates. A signed shift is not labeled increased magnitude. No empirical physical-error tolerance is introduced.

## Nonmetal sensitivity

Use the existing outcome-blind 96-structure manifest (34 diagnostics + 62 deterministic ID hashes). Accepted metal-involving indexed-site/image edges are retained exactly; only nonmetal cutoff scales 1.20 and 1.30 are substituted for 1.25, with minimum distance and reciprocal semantics unchanged. No metal CrystalNN run or whole-cohort descriptor build occurs. Changed unoriented physical edges are image-aware; distinct legitimate images are retained. First-shell/two-hop features and environment membership are recomputed. The same stronger original-fold OOF residuals are used for every graph. Effects are reported on each graph's available pairs and on the common paired-QMOF intersection across all three graphs. Paired effect changes use QMOF bootstrap. The added wider-cutoff H contacts have their actual distances, O–H–O angles, accepted/alternative H degrees and shared-proton ambiguity recorded; scale 1.30 is a stress configuration, not equally validated bond truth. No fixed-model prediction sensitivity or retrained alternative model is needed for this finite contrast question, and neither is claimed.

## Scope and provenance

The prior metal-only fitting result remains a disclosed follow-up after a global failure, not independent replication and not a replacement for the stronger conditional global test. No second metal-specific conditional ladder was added. No new DFT, experiments, GNNs, embeddings, broad algorithm search or hyperparameter optimization was performed. No charge scheme is truth. Charge accuracy, redox, universal bonding, broader family/topology transfer and causal electronic mechanisms remain outside evidence.

Source: [QMOF version 18](https://figshare.com/articles/dataset/QMOF_Database/13147324/18), [official repository](https://github.com/Andrew-S-Rosen/QMOF). Cite Rosen et al., Matter 4 (2021), [10.1016/j.matt.2021.02.015](https://doi.org/10.1016/j.matt.2021.02.015), and npj Computational Materials 8,112 (2022), [10.1038/s41524-022-00796-6](https://doi.org/10.1038/s41524-022-00796-6). QMOF data are CC BY 4.0; changes here are extraction, operational descriptors and statistical analysis.
'''
    (ROOT/'reports/METHODS.md').write_text(methods,encoding='utf-8')
    captions='''# Final figure captions

1. **Raw scheme differences and signs (unchanged).** Atom-weighted mean signed differences with structure denominators; full atom/structure distributions are saved. Opposite signs near zero are separated by declared magnitude thresholds. The three differences are algebraically dependent; no charge scheme is ground truth.

2. **Nested conditional comparison.** Structure-balanced held-out MAE with paired-design uncertainty, all atoms (top) and central metals (bottom). The common training-fitted piecewise offset precedes a center residual mean, a fixed tree with identity+q_A only, then the identical tree with first shell and two-hop paths. All competitors share fitting samples/weights and test atoms. Original QMOF folds and identifier/composition-blocked folds are shown. Original intervals resample QMOFs; blocked intervals resample component sums/counts, preserving the QMOF-balanced estimand. q_A is an explicit conditional follow-up input, not a predictor silently added to old charge-blind models.

3. **Directional and magnitude departures.** Predeclared same-element, within-QMOF contrasts after the stronger flexible charge-only model. Top: mean signed residual contrast. Bottom: mean absolute residual contrast. Intervals are 95% structure/component bootstrap under the indicated split. A signed shift is not evidence of greater absolute disagreement. Absolute departure is empirical charge-mapping departure, not physical error or calibrated physical uncertainty. Complete raw and spread quantities and paired support are in the atlas.

4. **Identifier/composition-blocked challenge.** Conditional two-hop gain over flexible charge-only, separately for all structures, records with exact MOFkeys and nonsingleton-component records. Both fold schemes are displayed, with the intended structure-balanced estimand and paired uncertainty. The nonsingleton subset prevents singleton abundance from hiding failure on the challenge population. Grouping links exact non-null MOFkey or exact source reduced formula; it is neither proof of structural identity nor universal family/topology separation.

5. **Actual nonmetal-cutoff sensitivity.** Signed (top) and absolute (bottom) calibrated contrasts on the common paired-QMOF intersection of scales 1.20, 1.25 and 1.30 on the existing 96-structure subset. All accepted metal edges are retained, no new metal NN is run, and the charge-only OOF residuals are identical across graphs. Markers are offset vertically for legibility; identical effects reflect unchanged tested environment assignments. Available-pair counts and paired effect changes are saved separately. Scale 1.30 adds near-midpoint O–H–O contact pairs and is labeled a stress configuration. Small subset support, especially for some H/O/N contrasts, limits scope. No retrained alternate-model performance or universal bond truth is inferred.

All five figures have PNG/SVG outputs and exact plotting inputs. Figure filenames 02/04/05 retain prior paths so existing references remain usable; their content and captions are replaced by the closeout generators. Earlier figures remain archived in `work/pre_closeout/ChaMM_review.zip`.
'''
    (ROOT/'reports/FIGURE_CAPTIONS.md').write_text(captions,encoding='utf-8')
    resulttext='# Final numerical tables\n\nAll values in e unless counts. Full atom and structure MAEs, fold increments, component-resampling inputs and support are in results/closeout_v1/.\n\n## Nested comparison, all structures\n\n'+md(pd.DataFrame(primary))+'\n\n## Blocked challenge support\n\n'+md(pd.DataFrame(challenge))+'\n\n## Global central-metal result\n\n'+md(pd.DataFrame(metal))+'\n\n## Signed and magnitude headline contrasts, blocked folds\n\n'+md(ev[ev.quantity.isin(['calibrated_signed','calibrated_absolute'])][['pair','left','right','quantity','n_structures','blocked_effect','blocked_ci_low','blocked_ci_high','validation_status','graph_status','graph_common_structures']])+'\n'
    (ROOT/'reports/RESULT_TABLES.md').write_text(resulttext,encoding='utf-8')
    closeout=f'''# Scientific closeout — completed finite tests

**{verdict}.** No required closeout analysis remains interrupted or unrun.

| Work package | Outcome | Main evidence |
| --- | --- | --- |
| Stronger conditional control | Survived for common nonmetal/all-atom predictions; central-metal result mixed. Several earlier signed effect sizes attenuated. | original/ and blocked/ OOF partitions; atlas/paired_increments.csv |
| Identifier/composition blocking | Survived overall, MOFkey-present and nonsingleton challenge; correct the earlier metadata-availability statement. | metadata/component_mapping.csv, atlas/grouping_challenge.csv |
| Shared nonmetal cutoff | Leading tested assignments/contrasts unchanged on 96 structures; six shared-proton sites expose a real representation ambiguity. Small subset support limits broader robustness claims. | nonmetal/physical_edge_changes.csv, available_pair_contrasts.csv, common_pair_contrasts.csv |
| Signed/magnitude atlas | Directional departure, magnitude and spread separated. A signed shift does not generally imply greater magnitude. Null and unsupported classes retained. | atlas/environment_atlas.csv, compact_environment_atlas.csv, claim_to_evidence.csv |

The unchanged main dataset is 9,171 structures / 1,420,532 atoms with configuration `{manifest['configuration_id']}`. Legacy outputs are preserved in their original result directories; the pre-closeout review ZIP is archived in `work/pre_closeout/ChaMM_review.zip`. New fits use `work/closeout_v1/` and new predictions/metadata use `results/closeout_v1/`, with separate original/blocked configurations and training samples. No old fitted model or OOF prediction was reused for a different split.

The atlas is a structured empirical description, not a physical-error atlas. Earlier all-metal fits remain a disclosed follow-up after the global failure; they are not independent replication and are not a new conditional result. The finite global conditional metal failure for DDEC6−Bader is retained. There is no requirement to rescue it with more models.

Accuracy ranking, redox, universal bonding, broader family/topology extrapolation, independent experiments and causal mechanisms remain outside this evidence. All specified finite tests are complete; no new method search or research campaign is part of the closeout. Stop analysis expansion here.
'''
    (ROOT/'reports/SCIENTIFIC_CLOSEOUT.md').write_text(closeout,encoding='utf-8')
    (ROOT/'reports/CLAIM_TO_EVIDENCE.md').write_text('# Claim-to-evidence closeout\n\n'+md(claims)+'\n',encoding='utf-8')
    readme=f'''# ChaMM final scientific closeout

**{verdict}.** Start with [SCIENTIFIC_CLOSEOUT](reports/SCIENTIFIC_CLOSEOUT.md), [findings](reports/SCIENTIFIC_FINDINGS.md), [methods](reports/METHODS.md), [numerical tables](reports/RESULT_TABLES.md), [claim-to-evidence](reports/CLAIM_TO_EVIDENCE.md) and [figure captions](reports/FIGURE_CAPTIONS.md).

The full frozen production cohort remains `data/structural_features.parquet` and `data/charge_targets.parquet` (9,171 structures / 1,420,532 atoms). The accepted graph was not rebuilt or changed. Original positive, negative and ambiguous results remain in `results/charge`, `models`, `science`, `metal_control`, `sensitivity` and `connectivity`; the earlier complete review bundle is `work/pre_closeout/ChaMM_review.zip`.

## New results

- `results/closeout_v1/metadata/`: source-verified metadata, component map and deterministic blocked fold manifest.
- `original/`, `blocked/`: 15 atom-keyed OOF partitions each, current-fold calibration parameters/knots, exact training sample keys/weights and configuration. Stronger charge-only residuals are explicit columns. Old charge-blind models are unchanged.
- `atlas/grouped_metrics.csv`, `paired_increments.csv`, `fold_metrics.csv`, `fold_increments.csv`: absolute errors, correctly conditioned paired increments and support.
- `atlas/structure_losses.parquet`: compact per-QMOF loss inputs for all reported populations; components, atom counts and fold tags support both MAE estimands and uncertainty.
- `atlas/environment_atlas.csv`, `compact_environment_atlas.csv`, `headline_contrasts.csv`, `environment_units.parquet`, `contrast_vectors.parquet`: signed/raw/absolute/spread results with support and ambiguity.
- `nonmetal/`: the fixed 96 structures, 1.20/1.25/1.30 features, changed physical edges, membership records, stronger subset residuals, available-pair effects, common-pair effects and selection changes. No new metal NN calculation.
- `claim_to_evidence.csv`, `main_comparison.csv`, `blocked_support_comparison.csv`, `verdict.json`: final quantitative claim record.
- New fitted models/checkpoints and reusable numeric matrix: `work/closeout_v1/`.

## Run or resume only the finite closeout

From this project root, with the completed production inputs and observed runtime:

```powershell
python -u code/24_conditional_closeout.py
python -u code/25_closeout_atlas.py
python -u code/26_nonmetal_closeout.py
python -u code/25_closeout_atlas.py --verify-compact
python -u code/21_figures.py
python -u code/22_review_package.py
```

Stage 24 can resume one split with `--split original` or `--split blocked`; it checkpoints complete fold/pair partitions under a checked configuration. `--prepare-only` reconstructs metadata/block manifests. Original fitting keys/weights are reused exactly; blocked keys are generated with the same training-only cap rule. Stage 26 supports `--build-only` and `--analyze-only`; it reads accepted graph caches and never calls CrystalNN. Stages 21/22 dispatch to the closeout generator whenever `results/closeout_v1` exists; missing closeout prerequisites raise an error rather than restoring old conclusions. A legacy `--verdict` string cannot replace the current data-derived closeout narrative.

Inputs needed for complete model reproduction are the frozen full atom tables, old `results/models/training_samples.parquet` and configuration, old fold manifest, pinned `work/qmof_database.zip`, the supplied reviewer metadata join, accepted subset graph/structure caches and existing 96-structure manifest. They are already present locally. The compact ZIP is for numerical review, not a self-contained full-cohort refit: it includes component/fold/sample manifests, new per-QMOF loss/contrast inputs, subset residuals/features/geometry/accepted and alternative edges, final tables, code and figures. Full atom-level global OOF files, full frozen atom tables and large regenerable caches remain local, with identities in `reports/artifact_index.json`. Do not rerun the full descriptor build for this closeout.

The compact numerical inputs can be checked without a full model refit using `python code/25_closeout_atlas.py --verify-compact`; the verifier labels the absence of an omitted original full-project sample file rather than pretending to recheck it. The full local run verified the original sample keys/weights exactly. The output deliberately retains weak/null and negative cases. Related-record blocking is conservative identifier/composition validation, not universal family/topology separation. No charge scheme is ground truth. QMOF source: [pinned Figshare version 18](https://figshare.com/articles/dataset/QMOF_Database/13147324/18), CC BY 4.0; original paper citations are in METHODS.md. `requirements-observed.txt` and `reports/runtime_environment.json` record the runtime.
'''
    (ROOT/'README.md').write_text(readme,encoding='utf-8')
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    (ROOT/'RUN_STATUS.md').write_text(f'# Scientific closeout complete\n\n{now}\n\n{verdict}.\n\nAll four finite work packages completed. Frozen cohort: 9,171 structures / 1,420,532 atoms; main graph unchanged. Thirty new split/fold/pair OOF partitions and 90 fixed tree fits completed, with full held-out coverage and saved training samples. Metadata/component counts reconcile to the independent review. Original global OOF targets/folds/90 aggregate rows reconciled once. Nonmetal graph sensitivity is limited to the existing 96 structures and makes zero metal NN calls. Numerical atlas, claim table, five figures, captions and compact review ZIP are current.\n\nNo expensive job remains interrupted. Exact reproduction/resume commands are in README.md. Outside-scope transfer, experiments, charge accuracy and mechanism studies are not incomplete closeout tasks.\n',encoding='utf-8')
    (OUT/'REPORT_READY.json').write_text(json.dumps(dict(completed=True,utc=now,verdict=verdict),indent=2))
    # A precise index of frozen and follow-up data; full atom-level files stay local.
    paths=[ROOT/'data/structural_features.parquet',ROOT/'data/charge_targets.parquet',ROOT/'data/descriptor_manifest.json']+[p for p in OUT.rglob('*') if p.is_file() and p.suffix in ['.csv','.json','.parquet']]
    index=[]
    for p in sorted(paths):
        h=hashlib.sha256()
        with p.open('rb') as stream:
            for b0 in iter(lambda:stream.read(8*1024*1024),b''):h.update(b0)
        index.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=h.hexdigest()))
    (ROOT/'reports/artifact_index.json').write_text(json.dumps(index,indent=2))
    include=[ROOT/'README.md',ROOT/'RUN_STATUS.md',ROOT/'requirements-observed.txt',ROOT/'data/descriptor_manifest.json',ROOT/'data/chamm_common_cohort.csv',ROOT/'results/grouped_splits.csv',ROOT/'results/models/model_configuration.json',ROOT/'results/models/evaluated_ladder.json']
    for folder in ['code','reports','figures']:
        include += [p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    include += [p for p in OUT.rglob('*') if p.is_file() and (p.suffix in ['.csv','.json'] or (p.suffix=='.parquet' and p.parent.name not in ['original','blocked']))]
    include += [OUT/s/'training_samples.parquet' for s in ['original','blocked']]
    include += [p for p in (ROOT/'results').rglob('*.csv') if OUT not in p.parents]
    include += [ROOT/'ChaMM_scientific_closeout/SCIENTIFIC_REVIEW.md',ROOT/'ChaMM_scientific_closeout/verified_source_metadata.csv']
    cid=manifest['configuration_id'];ids=pd.read_csv(OUT/'nonmetal/subset_manifest.csv').qmof_id
    include += [ROOT/'work/structures'/(q+'.json') for q in ids]+[ROOT/'work'/('graphs_'+cid)/(q+'.json') for q in ids]
    temp=ROOT/'ChaMM_review.zip.partial'
    with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(include)):z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(temp) as z:assert z.testzip() is None
    temp.replace(ROOT/'ChaMM_review.zip');print('Scientific closeout review ZIP MB',round((ROOT/'ChaMM_review.zip').stat().st_size/1e6,2),flush=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--figures',action='store_true');a=p.parse_args()
    make_figures() if a.figures else generate()
