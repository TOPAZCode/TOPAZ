all: mstw mrst2001lo cteq6 cteq10 convert

mrst: mrst2001lo.o
	ifort -O2 -c mrst2001lo.f

mstw: mstwpdf.o
	ifort -O2 -c mstwpdf.f

cteq6: Cteq66Pdf.o
	ifort -O2 -c Cteq66Pdf.f

cteq10: CT10Pdf.o
	ifort -O2 -c CT10Pdf.f

convert: cteq2mrst.o
	ifort -O2 -c cteq2mrst.f
