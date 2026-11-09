"""Check copied estimates, figure inputs, bibliography, and TeX dependencies.

This is a compact-package check, not a full refit or independent target audit.
"""
from pathlib import Path
import hashlib
import json
import re
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[2]

def main():
    checks=[]
    def check(condition,message):
        if not bool(condition): raise AssertionError(message)
        checks.append(message)
    inp=ROOT/'figures/inputs';ov=ROOT/'paper'
    for item in json.loads((inp/'PROVENANCE.json').read_text()):
        digest=hashlib.sha256((inp/item['file']).read_bytes()).hexdigest()
        check(digest==item['sha256'],'Frozen input matches source: '+item['file'])
    metrics=pd.read_csv(inp/'conditional_metrics.csv')
    gains=pd.read_csv(inp/'conditional_gains.csv')
    main=pd.read_csv(ROOT/'evidence/main_comparison.csv')
    for _,r in main.iterrows():
        for model in ['piecewise','center','charge_only','conditional_first','conditional_twohop']:
            row=metrics[metrics.split.eq(r['split'])&metrics.pair.eq(r.pair)&metrics.stratum.eq('all')&metrics.model.eq(model)].iloc[0]
            check(np.isclose(row.mae_structure,r[model],atol=1e-12,rtol=0),f'Main table reconciles: {r["split"]}/{r.pair}/{model}')
    metal=gains[gains.split.eq('blocked')&gains.pair.eq('ddec6_bader')&gains.stratum.eq('metal')&gains.comparison.eq('structure_given_q_A')].iloc[0]
    check(metal.ci_high<0,'Negative global conditional metal result retained')
    edge=pd.read_csv(inp/'nonmetal_edges.csv');change=pd.read_csv(inp/'nonmetal_changes.csv')
    check(len(edge)==12,'Twelve wider-cutoff physical edges')
    c=change[change.scale.eq(1.3)]
    check(c.first_shell_changed.sum()==18 and c.twohop_changed.sum()==38,'Nonmetal stress feature counts reconcile')
    bib=(ov/'references.bib').read_text();keys=set(re.findall(r'@\w+\s*\{\s*([^,]+)',bib))
    for f in [ov/'main.tex',ov/'si.tex',*sorted((ov/'sections').glob('*.tex'))]:
        s=f.read_text()
        for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}',s):
            for key in group.split(','):check(key.strip() in keys,f'Resolved citation {key.strip()} in {f.name}')
        for name in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',s):
            check((ov/name).exists(),f'Figure present: {name}')
        for name in re.findall(r'\\input\{([^}]+)\}',s):
            check((ov/(name+'.tex')).exists(),f'Section present: {name}')
    for fig in (ov/'figures').glob('*.pdf'):
        check(fig.with_suffix('.svg').exists() and fig.with_suffix('.png').exists(),f'Three formats: {fig.stem}')
    result={'status':'PASS','n_checks':len(checks),'checks':checks,
            'scope':'Exact source copies and compact estimate reconciliation; no full model refit.'}
    (ROOT/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS: {len(checks)} compact-package checks')

if __name__=='__main__':main()
