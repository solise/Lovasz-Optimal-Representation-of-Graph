function lanzador_PSOro_Int(fichero)
%% Lanzador del algoritmo PSO para el problema de Rango Ortogonal de un grafo

%Ejemplo para el Grafo C5

%n = 13; %Orden del grafo
%e = 39; %Tamaño del grafo
n = 8;e = 16;
theta = 8/(2+sqrt(2)); %Theta de lovasz
%Matriz de adyacencia
% M = [0 1 0 0 1;
%  1 0 1 0 0;
%  0 1 0 1 0;
%  0 0 1 0 1;
%  1 0 0 1 0];
%M=[0 1 1 0 0 0 1 1 0 0 0 1 1;1 0 1 1 0 0 0 1 1 0 0 0 1;1 1 0 1 1 0 0 0 1 1 0 0 0;0 1 1 0 1 1 0 0 0 1 1 0 0;0 0 1 1 0 1 1 0 0 0 1 1 0;0 0 0 1 1 0 1 1 0 0 0 1 1;1 0 0 0 1 1 0 1 1 0 0 0 1;1 1 0 0 0 1 1 0 1 1 0 0 0;0 1 1 0 0 0 1 1 0 1 1 0 0;0 0 1 1 0 0 0 1 1 0 1 1 0;0 0 0 1 1 0 0 0 1 1 0 1 1;1 0 0 0 1 1 0 0 0 1 1 0 1;1 1 0 0 0 1 1 0 0 0 1 1 0];
% Circulante C8
%M = [0 1 1 0 0 0 1 1;1 0 1 1 0 0 0 1; 1 1 0 1 1 0 0 0; 0 1 1 0 1 1 0 0; 0 0 1 1 0 1 1 0; 0 0 0 1 1 0 1 1; 1 0 0 0 1 1 0 1; 1 1 0 0 0 1 1 0]
M = load(fichero)
PSOroInt(n, M, theta, e)
end