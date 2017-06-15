function [phimenor,Sol,dimension] = PSOroFiel_valores(w, xi, n, M, e)
%Semejante a la versión PSOroFiel pero con parámetros de entrada dados (w,xi)

% PARÁMETROS DE ENTRADA:
%  - M = MATRIZ DE ADYACENCIA;
%  - n = Orden del grafo;
%  - e = Nº Aristas del grafo;
%  - theta = theta de Lovasz
%
% PSEUDOCÓDIGO:
% - Definir Funcion objetivo f(x).
% - Establecer restricciones del problema, lineales (lincon)
% y no lineales (nonlcon).
% - Establecer condiciones a aplicar sobre los puntos fuera
% de la region factible (outfac).
% - Determinar la dimension del enjambre (N) y numero de
% iteraciones (iter).
% - Inicializar posiciones y velocidades iniciales.
% - Para un numero 'iter' de iteraciones:
% - Para cada elemento i en N:
% * Chequear las restricciones: lincon(xi) y nonlcon(xi).
% * Si no cumple alguna de las restricciones
% -> Aplicar outfac(xi).
% * Calcular f(xi).
% * Actualizar mejor posicion, pi.
% * Actualizar velocidad, vi.
% * Actualizar posicion actual, xi.
% - Fin para.
% - Actualizar mejor posicion global, pg.
% - Fin para.

%%% Constantes %%%
iter = 500; %Número de iteraciones
N = 50; %Tamaño del enjambre
tolerancia=10e-6;

%w=0.4; %Factor de inercia
%xi=1.4; %Factor de constricción

%%%%%%%%%%%%%%%%%%%%%
%Número de dimensiones con las que probar
ncores=n-3;

GeneralPopulation = cell(ncores);
GeneralCosteMejor = cell(ncores);
GeneralphiMejor = cell(ncores);
GeneralMejor = cell(ncores);
Generald = cell(ncores); %dimensión

for j=1:ncores
    GeneralPopulation{j} = CalculaPoblacionInic(N,n,M,j+2);
    GeneralCosteMejor{j}=Inf;
    GeneralphiMejor{j} =Inf;
    GeneralMejor{j} = bee(n,j+2);
    Generald{j}=j+2;
end

%% Bucle principal

% INTENTA CONSTRUIR UNA REPRESENTACIÓN POR CADA DIMENSIÓN
for idim = 1:ncores
    
    % Toma valores locales para cada dimensión
    Population = GeneralPopulation{idim};
    phiMejor = GeneralphiMejor{idim};
    Mejor = GeneralMejor{idim};
    d = Generald{idim};
    
    MejorPop = Population;       %Lista con los mejores de cada generación de cada abeja
    PhiMejorPop = ones(1,N)*Inf; %Valores de mejor distancia obtenida por abeja
    IterSinMejora = 0;

    %Iteración de optimización 
    while phiMejor >= tolerancia && IterSinMejora <= iter       
        for j =1:N %Para cada abeja
            beej = Population(j);
            
            phij = f(beej.data(),M,d,n,e); %Calcula la función de coste sobre cada individuo según las condiciones de factible

            %Actualiza Valores locales de la abeja
            if phij < PhiMejorPop(j)
                MejorPop(j) = beej;
                PhiMejorPop(j) = phij;
            end
            
            %Actualiza valores globales para esta dimensión
            if phij < phiMejor
                if phiMejor - phij >= 10e-2
                     IterSinMejora = 0;
                end
                Mejor.setData(beej.data());
                phiMejor = phij;
            end

            %Movimiento de las abejas
            mbee = MejorPop(j);
            Vbee = beej.dataVel(); %Velocidad de la abeja
            r = unifrnd(0,1);
            Velocidad = w*Vbee + xi*(r*(Mejor.data() - beej.data()) + (1-r)*(mbee.data() - beej.data())); %Calcula el vector de velocidad
            beej.setVel(Velocidad);
            newDato = beej.data()+Velocidad;
            beej.setData(newDato/norm(newDato));
            Population(j) = beej; %Obtiene nueva situación del individuo

        end  
        
        if w > 0
            w = w-0.1; %Se reduce el factor de inercia
        end
        
        IterSinMejora = IterSinMejora +1;
    end
    
    %Almacenamiento en General
    GeneralPopulation{idim} = Population;
    GeneralphiMejor{idim} = phiMejor;
    GeneralMejor{idim} = Mejor;

end
% FIN DEL ANÁLISIS POR CADA DIMENSION


phimenor = Inf;
MejorDimension = 0;
MejorIndice=1;
for j=1:ncores
    if phimenor > GeneralphiMejor{j}
        phimenor = GeneralphiMejor{j};
        MejorDimension = j+2;
        MejorIndice=j;
    end
end

ElMejor = GeneralMejor{MejorIndice};


if phimenor < 10e-6
    disp('Encontrada Representación Óptima')
    dimension = MejorDimension
    Sol=ElMejor.data()
    phimenor
else
    disp('Mejor solución encontrada')
    dimension = MejorDimension
    Sol=ElMejor.data()
    phimenor
end

end