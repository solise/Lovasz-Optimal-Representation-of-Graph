$Ontext
   Resolución del problema de Representación Lovasz-optima

$Offtext

Scalar dim /9/;
Scalar vertices /20/;

Set
         i indice vertices /1*20/
         j indice dim /1*9/
         Alias(i,k)
         Alias(j,r)
;

Parameters

TABLE
 A(i,k)  Matriz de adyacencia de grafo

   1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
2  1  0  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
3  1  1  0  1  1  1  1  1  1  1  0  0  0  0  1  1  1  1  0  0
4  1  1  1  0  1  1  1  1  1  1  0  0  0  0  1  1  1  1  0  0
5  1  1  1  1  0  1  1  1  0  0  1  1  0  0  1  1  0  0  1  1
6  1  1  1  1  1  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1
7  1  1  1  1  1  1  0  1  0  0  0  0  1  1  0  0  1  1  1  1
8  1  1  1  1  1  1  1  0  0  0  0  0  1  1  0  0  1  1  1  1
9  1  1  1  1  0  0  0  0  0  1  1  1  1  1  1  1  1  1  0  0
10 1  1  1  1  0  0  0  0  1  0  1  1  1  1  1  1  1  1  0  0
11 1  1  0  0  1  1  0  0  1  1  0  1  1  1  1  1  0  0  1  1
12 1  1  0  0  1  1  0  0  1  1  1  0  1  1  1  1  0  0  1  1
13 1  1  0  0  0  0  1  1  1  1  1  1  0  1  0  0  1  1  1  1
14 1  1  0  0  0  0  1  1  1  1  1  1  1  0  0  0  1  1  1  1
15 0  0  1  1  1  1  0  0  1  1  1  1  0  0  0  1  1  1  1  1
16 0  0  1  1  1  1  0  0  1  1  1  1  0  0  1  0  1  1  1  1
17 0  0  1  1  0  0  1  1  1  1  0  0  1  1  1  1  0  1  1  1
18 0  0  1  1  0  0  1  1  1  1  0  0  1  1  1  1  1  0  1  1
19 0  0  0  0  1  1  1  1  0  0  1  1  1  1  1  1  1  1  0  1
20 0  0  0  0  1  1  1  1  0  0  1  1  1  1  1  1  1  1  1  0

;

Variables
x(i,j) Valor de la componente j-esima del vector i-esimo
y(j)   Handle
z      Valor de la Theta de Lovasz
;

*INTEGER VARIABLE x(i,j);
*INTEGER VARIABLE y(j);
x.up(i,j) = 1;
x.lo(i,j) = -1;
y.up(j) = 1;
y.lo(j) = -1;



Equations
    OBJ     funcion objetivo
    R1(i,k)   restriccion adyacencias
    R2(i)   restriccion norma 1
    R3      restriccion Handle
;

** La funcion objetivo es una restriccion de igualdad (=E=)

OBJ..            z =e= SUM(i,SUM(j,x(i,j)*y(j))*SUM(j,x(i,j)*y(j)));
R1(i,k)..        A(i,k)*SUM(j,x(i,j)*x(k,j)) =e= 0;
R2(i)..          SUM(j, x(i,j)*x(i,j)) =e= 1;
R3..             SUM(j, y(j)*y(j)) =e= 1;


** Definicion del modelo incluyendo todas las restricciones  declaradas
MODEL rolo /ALL/;
*schedule.optfile =1;
*** Opcion que indica el numero de iteraciones maxima
Option iterlim=1e9;
** Opcion que indica el tiempo maximo de ejecucion
*Option reslim=1e10;
** Opcion que indica que el optimizador se detiene cuando encuentra la solucion optima de manera exacta
Option optca=0;
Option optcr = 0;
option decimals = 8;

loop(i, loop(j, x.l(i,j) = 1 ) ) ;
loop(j, y.l(j)=1);

** solicita a GAMS que resuelva el problema mediante un optimizador
SOLVE rolo using nlp maximizing z;
** Muestra los valores de las variables
display  x.l, y.l, z.l;
