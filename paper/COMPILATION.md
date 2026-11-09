# Compilation

This directory is self-contained. In Overleaf, upload its contents and select
`main.tex`; select `si.tex` to compile the supporting information separately.
Use pdfLaTeX and BibTeX. Required figure PDF masters are in `figures/`.

Locally, from the repository root:

```sh
python reproduce/build_paper.py
```

Alternatively, in this directory run `pdflatex main`, `bibtex main`, and
`pdflatex main` twice; repeat with `si`. Both `.bbl` files are supplied for
portable compilation. Figure SVGs and PNGs are editable/preview exports;
the manuscript embeds only vector PDFs. This is a journal-neutral article
layout. Author metadata is preserved from the supplied manuscript.
