"""Generate main and SI numerical tables from saved estimates; no hand-entered values."""
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2]
D=pd.read_csv(R/'figures/inputs/conditional_metrics.csv')
G=pd.read_csv(R/'figures/inputs/conditional_gains.csv')
PAIRS=['ddec6_cm5','ddec6_bader','cm5_bader']
def name(p):return p.replace('ddec6','DDEC6').replace('cm5','CM5').replace('bader','Bader').replace('_','--')
def table(pop,caption,label):
 rows=[]
 for p in PAIRS:
  d=D[D.split.eq('blocked')&D.pair.eq(p)&D.stratum.eq(pop)].set_index('model')
  g=G[G.split.eq('blocked')&G.pair.eq(p)&G.stratum.eq(pop)&G.comparison.eq('structure_given_q_A')].iloc[0]
  vals=[d.loc[m,'mae_structure'] for m in ['charge_only','conditional_first','conditional_twohop']]
  rows.append(name(p)+' & '+' & '.join(f'{v:.5f}' for v in vals)+f' & {g.improvement_structure:.5f} & [{g.ci_low:.5f}, {g.ci_high:.5f}] '+r'\\')
 return '\n'.join([r'\begin{table}[!htbp]',r'\centering\small',r'\caption{'+caption+'}',r'\label{'+label+'}',r'\setlength{\tabcolsep}{3pt}',r'\begin{tabular}{lrrrrl}',r'\toprule',r'Charge pair & Charge only & First shell & Two-hop & Gain & 95\% interval \\',r'\midrule',*rows,r'\bottomrule',r'\end{tabular}',r'\end{table}',''])
(R/'paper/sections/main_table.tex').write_text(table('all',r'Identifier/composition-blocked structure-balanced MAE and paired two-hop gain over charge only, in e. All 9,171 structures and 1,420,532 atoms are evaluated. Intervals are pointwise 95\% component-bootstrap intervals conditional on fitted models.','tab:main'))
(R/'paper/sections/si_metal_table.tex').write_text(table('metal',r'Central-metal MAE and paired two-hop gain under identifier/composition blocking, in e. The global models are evaluated on 46,994 metal atoms in 9,165 structures. Pointwise 95\% component-bootstrap intervals condition on the fitted models.','tab:si-metal'))
print('Regenerated main all-atom and SI central-metal numerical tables.')
