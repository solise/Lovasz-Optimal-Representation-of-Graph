function lanzador_PSOro_valores(fichero,n,e,theta)
%% Lanzador del algoritmo PSO para el problema de Rango Ortogonal de un grafo

%Ejemplo para el Grafo C5
format longE
%n = 13; %Orden del grafo
%e = 39; %Tamaño del grafo
%% CIRCULANTE C8(1,2)
% n = 8;e = 16;
% theta = 8/(2+sqrt(2)); %Theta de lovasz
%% GRAFO CS3
% n=10;e=24;
% theta = 3.0945;
%% Matriz de adyacencia
% Circulante C8(1,2)
%M = [0 1 1 0 0 0 1 1;1 0 1 1 0 0 0 1; 1 1 0 1 1 0 0 0; 0 1 1 0 1 1 0 0; 0 0 1 1 0 1 1 0; 0 0 0 1 1 0 1 1; 1 0 0 0 1 1 0 1; 1 1 0 0 0 1 1 0]

M=load(fichero)

Costestxt = strcat('Tabla_Costes_500_',fichero,'.txt');
Mejorestxt = strcat('Mejores_PSO_500_',fichero,'.txt');

w=0.1;
while w <= 1
    xi = 1;
    while xi <= 3
        [CosteMejor, Mejor,dimension] = PSOroParallel_valores(w,xi,n, M, theta, e);
        Tupla = [w, xi, CosteMejor, dimension];
        save(Costestxt,'Tupla','-ascii','-double','-append');
        if CosteMejor < 1*10e-2
            save(Mejorestxt,'CosteMejor','-ascii','-double','-append');
            save(Mejorestxt,'dimension','-ascii','-double','-append');
            save(Mejorestxt,'Mejor','-ascii','-double','-append');
        end
        xi = xi+0.1;
    end
    w = w+0.1;
end

end