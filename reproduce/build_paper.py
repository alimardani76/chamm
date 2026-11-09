"""Compile both documents with pdfLaTeX/BibTeX and reject unresolved references."""
from pathlib import Path
import subprocess,re
R=Path(__file__).resolve().parents[1]
for doc in ['main','si']:
 for cmd in [['pdflatex','-interaction=nonstopmode','-halt-on-error',doc+'.tex'],['bibtex',doc],['pdflatex','-interaction=nonstopmode','-halt-on-error',doc+'.tex'],['pdflatex','-interaction=nonstopmode','-halt-on-error',doc+'.tex']]:
  p=subprocess.run(cmd,cwd=R/'paper',stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  if p.returncode:raise RuntimeError(p.stdout)
 log=(R/'paper'/f'{doc}.log').read_text(errors='replace')
 assert not re.search(r'undefined|Overfull \\hbox|Overfull \\vbox',log,re.I),doc+' has unresolved references or overflow'
 print('PASS:',doc+'.pdf compiled with resolved citations/references and no overflow')
