OBJ=lb-memparcmp
all: ${OBJ}.pdf

${OBJ}.pdf: ${OBJ}.tex
	pdflatex ${OBJ}
	bibtex ${OBJ}
	pdflatex ${OBJ}
	pdflatex ${OBJ}

