Manuscript LaTeX skeleton for npj Microgravity

Files added:
- npj_microgravity.tex        (main LaTeX driver)
- manuscript.bib              (bibliography placeholder)
- sections/*.tex              (abstract, introduction, methods, results, discussion, acknowledgements)

How to build
- Recommended: install a TeX distribution (TeX Live or MikTeX).
- From the manuscript/ folder run (example):
  pdflatex npj_microgravity.tex
  bibtex npj_microgravity
  pdflatex npj_microgravity.tex
  pdflatex npj_microgravity.tex

Or use latexmk for convenience:
  latexmk -pdf npj_microgravity.tex

Notes and next steps
- Replace author names, affiliations, and section text with the final manuscript content.
- Place figure files in astronaut-opposite-forcing/figures/ with names such as figure1.png, figure2.png or update the paths in the main tex file.
- If npj provides an official LaTeX template (class/style files), replace the documentclass line in npj_microgravity.tex with the journal-provided class and adjust bibliography style accordingly.
