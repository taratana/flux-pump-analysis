reset
set terminal pdfcairo font "Times-New-Roman,18" enhanced 
set output "i2.pdf"

set datafile separator ","

set grid
set key right bottom
set key opaque box
set xlabel "{/Times-New-Roman-Italic=25 t} [s]"
set ylabel "{/Times-New-Roman-Italic=25 i}_{2} [A]" 
set xrange [0:1000]


plot "with-filter-i2-120.csv" using 1:2  with lines lw 3   title "120",\
	 "with-filter-i2-90.csv" using 1:2  with lines lw 3   title "90",\
	 "with-filter-i2-60.csv" using 1:2  with lines lw 3   title "60",\
	 "with-filter-i2-30.csv" using 1:2  with lines lw 3   title "30",\
	 "with-filter-i2-0.csv" using 1:2  with lines lw 3   title "0",\


############################################################

