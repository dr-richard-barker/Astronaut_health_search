# LaTeX manuscript (npj Microgravity / Springer Nature style)

Assembles the opposite-forcing manuscript in the official **Springer Nature LaTeX
template** (`sn-jnl` class) — the format npj Microgravity accepts and typesets
from.

```
latex/
├── main.tex          # assembled manuscript (sn-jnl class, sn-nature refs)
├── references.bib    # 10 references transcribed from the manuscript
├── figures/          # the 13 figures (F1–F12 + F4b, PNG)
└── README.md         # this file
```

## How to compile

`sn-jnl.cls` / `sn-nature.bst` ship with Springer Nature's official template
(not vendored here).

- **Overleaf (recommended):** new project from the **"Springer Nature Article
  Template (sn-jnl)"** → replace its `main.tex` with this one, upload
  `references.bib` and `figures/`, compile with **pdfLaTeX**.
- **Local:** place `sn-jnl.cls` + `sn-nature.bst` here, then
  `pdflatex main` → `bibtex main` → `pdflatex main` → `pdflatex main`.

## Status / TODO before submission

- [ ] **Not yet compile-tested** — authored without a local TeX install; build
      once on Overleaf and fix any stragglers.
- [ ] **Author block** — co-authors and affiliation are `TODO` in `main.tex`
      (the source marks them "co-authors TBD", affiliation TODO).
- [ ] **Tables 1–13** — the manuscript references numbered tables (T1–T13, incl.
      sub-tables). The machine-readable CSVs are in `../../tables/`; typeset the
      key ones as LaTeX tables and/or bundle the rest as supplementary data per
      journal instructions. (Only figures are inlined here.)
- [ ] **References** — transcribed with volumes/pages; add DOIs (one entry,
      Duan 2023, needs volume/article-number — marked `% TODO`).
- [ ] **Figures** — repo PNGs; SVG editable versions exist in `../../figures/`.
      npj prefers vector/≥300 dpi for final submission.

## Source

Ported from `../manuscript.md`. Body text, figures, and references are the
author's own content — nothing was invented.
