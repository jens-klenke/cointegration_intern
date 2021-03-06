---
title: 'P-Approximation'
author: 'Jens Klenke and Janine Langerbein'
subtitle: 'Seminar in Econometrics'
type: "Term Paper"
discipline: "VWL M.Sc."
date: "today"
studid: "3071594"
supervisor: "Christoph Hanck"
secondsupervisor: "NA"
semester: "5"
estdegree_emester: "Winter Term 2020"
deadline: "Jan. 17th 2020"
output:
  pdf_document:
    keep_tex: yes
    template: template.tex
    fig_caption: yes
    citation_package: biblatex
    number_sections: true
toc: true
lot: true
lof: true
graphics: true
biblio-title: References
linkcolor: black
urlcolor: black
citecolor: black
colorlinks: true
font: Times New Roman
fontsize: 12pt
geometry: lmargin = 5cm,rmargin = 2.5cm,tmargin = 2.5cm,bmargin = 2.5cm
biblio-files: references.bib
classoption: a4paper
---

<!-- % Template Version 1.1 -->

<!-- Git version -->



# Introduction
Meta tests have been shown to be a powerful tool when testing for the null of non-cointegration. The distribution of their test statistic, however, is not always available in closed form. 

Start of the paper

## int 2 

# main 


@tibshirani_regression_1996


<!-- end of main part -->

\pagebreak

\pagenumbering{Roman}

\addcontentsline{toc}{section}{References}
\printbibliography[title = References]
\cleardoublepage

\begin{refsection}
\nocite{R-base}
\nocite{R-stargazer}
\nocite{R-stringr}
\nocite{R-tidyr}
\nocite{R-dplyr}
\nocite{R-glmnet}
\nocite{R-class}
\nocite{R-MASS}
\nocite{R-plm}
\nocite{R-leaps}
\nocite{R-caret}
\nocite{R-tree}
\nocite{R-gbm}
\nocite{R-plotmo}
\nocite{R-pls}
\nocite{R-splines}
\nocite{R-tictoc}
\nocite{R-plotly}
\nocite{R-inspectdf}
\nocite{R-rpart}
\nocite{R-rpart.plot}
\nocite{R-stargazer}
\nocite{R-knitr}
\nocite{R-purrr}
\nocite{R-randomForest}
\nocite{R-rstudioapi}





\nocite{R-Studio}

\printbibliography[title = Software-References]
\addcontentsline{toc}{section}{Software-References}
\end{refsection}


<!---
--------------------------------------------------------------------------------
------------- Appendix ---------------------------------------------------------
--------------------------------------------------------------------------------
-->

\cleardoublepage
\appendix
\setcounter{table}{0}
\setcounter{figure}{0}
\renewcommand{\thetable}{A\arabic{table}}
\renewcommand{\thefigure}{A\arabic{figure}}

\newgeometry{top = 2.5cm, left = 5cm, right = 2.5cm, bottom = 2cm}

# Appendices

Better sorting of the Appendix 


<!-- 
--------------------------------------------------------------------------------
------------- End of Appendix --------------------------------------------------
--------------------------------------------------------------------------------
--> 
\restoregeometry

\cleardoublepage

