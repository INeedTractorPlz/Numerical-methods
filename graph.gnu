#! /usr/bin/gnuplot -persist
set terminal jpeg
set output "alpha.jpg"
plot 'file_alpha.dat' using 1:6 w l,\
'file_alpha.dat' using 1:7 w l,\
'file_alpha.dat' using 1:8 w l,\
'file_alpha.dat' using 1:9 w l,\
'file_alpha.dat' using 1:10 w l,\
'file_alpha.dat' using 1:11 w l,\
#'file_alpha.dat' using 1:12 w l,\
#'file_alpha.dat' using 1:13 w l,\
#'file_alpha.dat' using 1:14 w l,\
#'file_alpha.dat' using 1:15 w l,\
#'file_alpha.dat' using 1:16 w l

