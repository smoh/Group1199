LATEX       = pdflatex -interaction=nonstopmode -halt-on-error
BASH        = bash -c
ECHO        = echo
RM          = rm -rf
TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out brf
CHECK_RERUN =

NAME = ms

all: ${NAME}.pdf


${NAME}.pdf: ${NAME}.tex ref.bib kk.tex
	${LATEX} ${NAME}
	bibtex ${NAME}
	${LATEX} ${NAME}
	( grep "Rerun to get" ${NAME}.log && ${LATEX} ${NAME} ) || echo "Done."
	( grep "Rerun to get" ${NAME}.log && ${LATEX} ${NAME} ) || echo "Done."

ref: ref.local.bib
	${RM} ref.bib;
	bibexport -a -o ref.bib ref.local

count: ms.tex
	texcount ms.tex

clean:
	${RM} $(foreach suff, ${TMP_SUFFS}, ${NAME}.${suff})
	${RM} gitstuff.aux
	${RM} *Notes.bib

submit: ${NAME}.pdf
	make
	tar -czvf archive.tar.gz $$(git ls-tree --name-only HEAD | grep -v ".gitignore\|Makefile\|figures") -C ./figures/ $$(cd figures && git ls-files)
	echo "archive.tar.gz created"

#NOTE:Pandoc does not run bibtex automatically (See https://github.com/jgm/pandoc/issues/1193)
report.pdf : report.tex
	# pandoc --natbib --template response.latex --variable=geometry:margin=1in --variable=fontsize:12pt -s -o report.tex report.md
	${LATEX} report
	bibtex report
	${LATEX} report
	( grep "Rerun to get" report.log && ${LATEX} report ) || echo "Done."
	( grep "Rerun to get" report.log && ${LATEX} report ) || echo "Done."

diff.pdf: ms.pdf submission/ms.tex
	latexdiff --disable-citation-markup -t BOLD submission/ms.tex ms.tex > diff.tex
	${LATEX} diff
	bibtex diff
	${LATEX} diff
	( grep "Rerun to get" diff.log && ${LATEX} diff ) || echo "Done."
	( grep "Rerun to get" diff.log && ${LATEX} diff ) || echo "Done."
