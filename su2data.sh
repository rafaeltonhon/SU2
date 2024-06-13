#!/usr/bin/bash
path=runs/beta=2.5/l=20
a=7
b=8
nc=1
nr=20
nt=3

cp su2treatdata.py su2observables.py su2statistics.py $path #runs/beta=2.5/l=20]
cp su2wilson.p su2pot.p su2gluon.p $path
cd $path
echo $nc > 'su2treatdata-in.dat'
echo $nt >> 'su2treatdata-in.dat'
echo $nr >> 'su2treatdata-in.dat'
echo $a >> 'su2treatdata-in.dat'
echo $b >> 'su2treatdata-in.dat'

python3 su2treatdata.py < 'su2treatdata-in.dat' # function ncorr nt nr a b
gnuplot su2pot.p
gnuplot su2gluon.p
gnuplot su2wilson.p
for i in $(seq 1 $a)
do 
	rm uncfort.$((i+100))
	rm uncfort.$((i+200))
	rm uncfort.$((i+300))
	
	mv meanfort.$((i+100)) data
	mv meanfort.$((i+200)) data
	mv meanfort.$((i+300)) data
done

pdflatex gluonpropagator.tex
pdflatex potentials.tex
pdflatex wilson.tex
pdflatex creutz.tex
pdflatex gluonformfactor.tex
find . -name "*.eps" -type f -delete
find . -name "*.tex" -type f -delete
find . -name "*.log" -type f -delete
find . -name "*.aux" -type f -delete
rm su2treatdata.py su2observables.py su2statistics.py su2wilson.p su2pot.p su2gluon.p
