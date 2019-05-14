reset
set terminal pdfcairo font "Times-New-Roman,18" enhanced 
set output "i_theta-t_curve.pdf"

set datafile separator ","

set key right bottom
set xlabel "{/Times-New-Roman-Italic=25 t} [s]"
set ylabel "{/Times-New-Roman-Italic=25 i}_{theta} [A]" 
set xrange [0:2000]


plot "output-without-filter.csv" using 1:2  with lines lw 3   title "without filter",\
	 "output-with-filter.csv" using 1:2  with lines lw 3   title "with filter"