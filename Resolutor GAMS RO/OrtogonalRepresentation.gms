$Ontext
   Resolución del problema de Representación ortogonal NO FIEL

$Offtext

Scalar dim /4/;
Scalar vertices /18/;

Set
         i indice vertices /1*18/
         j indice dim /1*4/
         Alias(i,k)
         Alias(j,r)
;

Parameters


TABLE
 A(i,k)  Matriz de adyacencia de grafo
    1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18
1   0  1  1  0  1  1  1  0  1  1  0  0  0  0  0  0  0  0
2   1  0  1  1  0  1  0  1  0  1  0  1  0  0  0  0  0  0
3   1  1  0  1  1  0  0  1  1  0  1  0  0  0  0  0  0  0
4   0  1  1  0  0  0  1  1  0  0  0  0  0  1  1  0  0  1
5   1  0  1  0  0  0  0  0  1  0  0  1  0  0  1  0  1  1
6   1  1  0  0  0  0  0  0  0  1  1  0  0  1  0  1  0  1
7   1  0  0  1  0  0  0  1  0  0  0  0  0  1  1  1  1  0
8   0  1  1  1  0  0  1  0  0  0  0  0  1  0  0  1  1  0
9   1  0  1  0  1  0  0  0  0  0  0  1  1  1  0  1  0  0
10  1  1  0  0  0  1  0  0  0  0  1  0  1  0  1  0  1  0
11  0  0  1  0  0  1  0  0  0  1  0  0  1  0  1  1  0  1
12  0  1  0  0  1  0  0  0  1  0  0  0  1  1  0  0  1  1
13  0  0  0  0  0  0  0  1  1  1  1  1  0  1  1  0  0  0
14  0  0  0  1  0  1  1  0  1  0  0  1  1  0  1  0  0  0
15  0  0  0  1  1  0  1  0  0  1  1  0  1  1  0  0  0  0
16  0  0  0  0  0  1  1  1  1  0  1  0  0  0  0  0  1  1
17  0  0  0  0  1  0  1  1  0  1  0  1  0  0  0  1  0  1
18  0  0  0  1  1  1  0  0  0  0  1  1  0  0  0  1  1  0

;


Variables
x(i,j)   Valor de la componente j-esima del vector i-esimo
h1(i,k)  Valor del producto de los no adyacentes  (debe ser distinto de 0)
z
;

POSITIVE VARIABLE h1(i,k);
x.up(i,j) = 1;
x.lo(i,j) = -1;

z.l = 1;


Equations
    OBJ       objetivo
    R1(i,k)   restriccion adyacencias
    R2(i,k)   restriccion no adyacencias
    R3(i)     restriccion norma 1
;

** La funcion objetivo es una restriccion de igualdad (=E=)

OBJ..                               z =e= 1;
R1(i,k)..                           A(i,k)*SUM(j,x(i,j)*x(k,j)) =e= 0;
R2(i,k)$(ord(i) lt ord(k))..        h1(i,k) =e= abs(mod(A(i,k)+1,2)*SUM(j,x(i,j)*x(k,j))) + A(i,k) ;
*R2(i,k)..        abs(mod(A(i,k)+1,2)*SUM(j,x(i,j)*x(k,j))) =g= mod(A(i,k)+1,2)*0.00000001;
R3(i)..                             SUM(j, x(i,j)*x(i,j)) =e= 1;


** Definicion del modelo incluyendo todas las restricciones  declaradas
MODEL ro /ALL/;
*schedule.optfile =1;
** Opcion que indica el numero de iteraciones maxima
Option iterlim=1e9;
** Opcion que indica el tiempo maximo de ejecucion
Option reslim=1e10;
** Opcion que indica que el optimizador se detiene cuando encuentra la solucion optima de manera exacta
Option optca=0;
Option optcr = 0;
option decimals = 8;

loop(i, loop(j, x.l(i,j) = 1 ) ) ;
*loop(i, loop(k, h1.l(i,k)$(ord(i) ge ord(k)) = 1) );  (* NO ES NECESARIO *)

** solicita a GAMS que resuelva el problema mediante un optimizador
SOLVE ro using dnlp minimizing z;
** Muestra los valores de las variables
display  x.l, z.l, h1.l;
