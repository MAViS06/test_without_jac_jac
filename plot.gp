splot 'trajectoire.log' w l ;
pause -1 "Appuyer sur Return pour continuer"

plot 'total_error.log' w l ;
pause -1 "Appuyer sur Return pour continuer"

plot 'lengths.log' u 1:2 w l, 'lengths.log' u 1:3 w l, 'lengths.log' u 1:4 w l ;
pause -1 "Appuyer sur Return pour continuer"

plot 'trajectoire.log' u 1:2 w l ;
pause -1 "Appuyer sur Return pour continuer"

plot 'trajectoire.log' u 1:3 w l ;
pause -1 "Appuyer sur Return pour continuer"

plot 'trajectoire.log' u 2:3 w l
pause -1 "Appuyer sur Return pour continuer"

