function sol =lanzador_PSOro(fichero)
%% Lanzador del algoritmo PSO para el problema de Rango Ortogonal de un grafo

%Ejemplo para el Grafo C5
format longE
%n = 13; %Orden del grafo
%e = 39; %Tama�o del grafo
% n = 8;e = 16;
% theta = 8/(2+sqrt(2)); %Theta de lovasz
n=10;e=24;
theta=3.0945;
%Matriz de adyacencia
% M = [0 1 0 0 1;
%  1 0 1 0 0;
%  0 1 0 1 0;
%  0 0 1 0 1;
%  1 0 0 1 0];
%M=[0 1 1 0 0 0 1 1 0 0 0 1 1;1 0 1 1 0 0 0 1 1 0 0 0 1;1 1 0 1 1 0 0 0 1 1 0 0 0;0 1 1 0 1 1 0 0 0 1 1 0 0;0 0 1 1 0 1 1 0 0 0 1 1 0;0 0 0 1 1 0 1 1 0 0 0 1 1;1 0 0 0 1 1 0 1 1 0 0 0 1;1 1 0 0 0 1 1 0 1 1 0 0 0;0 1 1 0 0 0 1 1 0 1 1 0 0;0 0 1 1 0 0 0 1 1 0 1 1 0;0 0 0 1 1 0 0 0 1 1 0 1 1;1 0 0 0 1 1 0 0 0 1 1 0 1;1 1 0 0 0 1 1 0 0 0 1 1 0];
% Circulante C8
%M = [0 1 1 0 0 0 1 1;1 0 1 1 0 0 0 1; 1 1 0 1 1 0 0 0; 0 1 1 0 1 1 0 0; 0 0 1 1 0 1 1 0; 0 0 0 1 1 0 1 1; 1 0 0 0 1 1 0 1; 1 1 0 0 0 1 1 0]
M=load(fichero)

sol = PSOroParallel(n, M, theta, e)
end