(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("sfmath" "cm") ("caption" "labelfont=bf") ("titlesec" "compact") ("xcolor" "table") ("todonotes" "colorinlistoftodos" "textsize=tiny") ("geometry" "margin=1in") ("biblatex" "backend=bibtex" "natbib=true" "style=numeric-comp" "sorting=none")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "sfmath"
    "caption"
    "subcaption"
    "url"
    "bm"
    "textpos"
    "graphicx"
    "titlesec"
    "amsmath"
    "esint"
    "amssymb"
    "textcomp"
    "enumitem"
    "fancyhdr"
    "hyperref"
    "booktabs"
    "tabularx"
    "tocloft"
    "setspace"
    "longtable"
    "xcite"
    "xcolor"
    "todonotes"
    "geometry"
    "biblatex"
    "placeins"
    "wrapfig"
    "ifthen"
    "pdfpages")
   (TeX-add-symbols
    "alox"
    "C")
   (LaTeX-add-bibliographies
    "library")))

