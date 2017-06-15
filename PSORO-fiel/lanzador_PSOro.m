function sol =lanzador_PSOro(fichero)
%% Lanzador del algoritmo PSO para el problema de Rango Ortogonal de un grafo

%Ejemplo para el Grafo C5
format longE
%n : Orden del grafo
%e : Tamaño del grafo

M=load(fichero)

n=size(M,1);
e=sum(sum(M))/2;

sol = PSOroFiel(n, M, e)
end