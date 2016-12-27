
$set n 1000

OPTIONS

OptCA = 0.0001
iterlim = 500000
work = 10000000
reslim = 50000;

FILE results /output.txt/

SETS

    i the rows /1 * 7/
    j the columns /1 * 7/
    k the color in the grid /1 * 8/
    ALIAS(i,m)
    ALIAS(j,n)
VARIABLES z;

BINARY VARIABLES C(i,j,k), CI(i,j,k);

EQUATIONS obj, position(i,j),middle(i,j,k),end(i,j,k);
    
    obj.. z =e= 1;

    position(i,j).. SUM(k,C(i,j,k) + CI(i,j,k)) =e= 1;
    
    middle(i,j,k).. SUM((m,n)$(C(i,j,k) eq 1),C(m,n,k) + CI(m,n,k)) =e= 3;
    end(i,j,k).. SUM((m,n)$(CI(i,j,k) eq 1),C(m,n,k) + CI(m,n,k)) =e= 2;

*   middle(i,j,k).. SUM((m,n)$(ord(m) le ord(i) + 1 and ord(m) ge ord(i) - 1 and ord(n) le ord(j) + 1 and ord(n) ge ord(j) - 1 and C(i,j,k) eq 1),C(m,n,k) + CI(m,n,k)) =e= 3;
*   end(i,j,k).. SUM((m,n)$(ord(m) le ord(i) + 1 and ord(m) ge ord(i) - 1 and ord(n) le ord(j) + 1 and ord(n) ge ord(j) - 1 and CI(i,j,k) eq 1),C(m,n,k) + CI(m,n,k)) =e= 2;




    CI.fx("1","7","8") = 1;
    CI.fx("2","2","2") = 1;
    CI.fx("2","3","1") = 1;
    CI.fx("2","6","4") = 1;
    CI.fx("3","1","3") = 1;
    CI.fx("3","4","3") = 1;
    CI.fx("3","5","6") = 1;
    CI.fx("3","6","5") = 1;
    CI.fx("5","4","5") = 1;
    CI.fx("6","2","1") = 1;
    CI.fx("6","6","7") = 1;
    CI.fx("7","1","2") = 1;
    CI.fx("7","2","4") = 1;
    CI.fx("7","4","8") = 1;
    CI.fx("7","6","8") = 1;
    CI.fx("7","7","7") = 1;

MODEL dots /all/;

SOLVE dots USING MIP
    
    MINIMIZING z;

    PUT results;


    LOOP(i, 
    	LOOP(j,
	    LOOP(k,
	    	PUT $(C.l(i,j,k) + CI.l(i,j,k) eq 1) k.tl:0:0, ", ";
		)
	    ) PUT /;
	);
    PUT /;
    PUT "modelstat= " dots.modelstat;
