set terminal pngcairo
set view equal xyz # to get an orthographic projection with all of the axes scaled equally
set xyplane at 0 # adjusts the position at which the xy plane is drawn in a 3D plot
set output "t.png"
set key left top
set title "Parametric plane with internal rays of main cardioid angle = 1/3 and  for period 3 component angle = 0"
set xlabel "cx"
set ylabel "cy"
# points for ray
set style line 1 lc rgb 'black' pt 7  ps 0.3  # ls 1 : black circle
# Plotting single points
set style line 2 lc rgb 'red' pt 7   ps 1.0 # ls 2 : red circle = root
set style line 3 lc rgb 'blue' pt 7   ps 1.0 # ls 3 : blue circle = center 0
set style line 4 lc rgb 'green' pt 7   ps 1.0 # ls 4 : green circle = center 3
# parametric curves 
set parametric
set trange [0:2*pi]
# main cardioid
fx(t) = cos(t)/2-cos(2*t)/4
fy(t) = sin(t)/2-sin(2*t)/4
# period 2 component
gx(t) = cos(t)/4-1
gy(t) = sin(t)/4
#
plot fx(t),fy(t) title "cardioid",\
	gx(t), gy(t) title "period 2",\
	 "m.txt" using 1:2 with points ls 1,\
	 '-' w p ls 3 title "center 0", '-' w p ls 4 title "center 3", '-' w p ls 2 title "root"
0.0 0.0
e
-0.122561166876654 0.744861766619744
e
-0.1249999999999998 0.6495190528383290
e

