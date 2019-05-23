reset
set terminal pdfcairo font "Times-New-Roman,18" enhanced 
set output "f2.pdf"

set datafile separator ","

set grid
set key right bottom
set key opaque box
set xlabel "{/Times-New-Roman-Italic=25 t} [s]"
set ylabel "{/Times-New-Roman-Italic=25 i}_{2} [A]" 
set xrange [0:1000]


plot "with-filter-f2-100.csv" using 1:2  with lines lw 3   title "{/=14 100 Hz}",\
	 "with-filter-f2-50.csv" using 1:2  with lines lw 3   title "{/=14 50 Hz}",\
	 "with-filter-f2-10.csv" using 1:2  with lines lw 3   title "{/=14 10 Hz}",\
	 "with-filter-f2-2.0.csv" using 1:2  with lines lw 3   title "{/=14 2 Hz}",\
	 "with-filter-f2-0.5.csv" using 1:2  with lines lw 3   title "{/=14 0.5 Hz}",\


############################################################

