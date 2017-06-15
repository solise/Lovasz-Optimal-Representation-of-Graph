function [ElCosteMejor, best,dimension] = PSOroParallel_valores(w, xi, n, M, theta, e)
% PAR�METROS DE ENTRADA:
%  - M = MATRIZ DE ADYACENCIA;
%  - n = Orden del grafo;
%  - e = N� Aristas del grafo;
%  - theta = theta de Lovasz
%
% PSEUDOC�DIGO:
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
iter = 500; %N�mero de iteraciones
N = 50; %Tama�o del enjambre
tolerancia=10e-15;

%w=0.4; %Factor de inercia
%xi=1.4; %Factor de constricci�n

%%%%%%%%%%%%%%%%%%%%%
%Obtiene el n�mero de nucleos
%ncores = feature('numCores');
ncores=n-3;

%Crea una piscina de trabajadores seg�n el n�mero de n�cleos
matlabpool(ncores);

Population = Composite(ncores);
CosteMejor = Composite(ncores);
phiMejor = Composite(ncores);
Mejor = Composite(ncores);
d = Composite(ncores); %dimensi�n

for j=1:ncores
    Population{j} = CalculaPoblacionInic(N,n,M,j+2);
    CosteMejor{j}=Inf;
    phiMejor{j} =Inf;
    Mejor{j} = bee(n,j+2);
    d{j}=j+2;
end

%% Bucle principal
%while (d<n) && (CosteMejor > 10e-6)
tic
spmd(ncores) 

    %Population = CalculaPoblacionInic(n,M,d);
    MejorPop = Population; %Lista con los mejores de cada generaci�n de cada abeja
    CosteMejorPop = ones(1,N)*Inf; %Costes de la mejor poblacion
    %CosteMejor=Inf; %Coste del mejor individuo
    PhiMejorPop = ones(1,N)*Inf; %Valores de mejor distancia obtenida por individuo
    %CosteMejor=Inf; %Coste del mejor individuo
    %phiMejor =Inf;
    %Mejor = bee(n,d);
    IterSinMejora = 0;

    %for i=1:iter %Iteraciones sobre cada dimensi�n
    while CosteMejor > 10e-2 && IterSinMejora <= iter        
        for j =1:N %Para cada abeja
            beej = Population(j);

            phij = f(beej.data(),M,d,n,e); %Calcula la funci�n de coste sobre cada individuo seg�n las condiciones de factible

            %Actualiza Valores locales de la abeja
            if phij < PhiMejorPop(j)
                MejorPop(j) = beej;
                PhiMejorPop(j) = phij;
            end
            if phij<=tolerancia && phiMejor <= tolerancia
                 Scorej = cost(beej.data(),theta,n,d); %Calcula la funci�n coste de optimizaci�n principal
                 if Scorej <= CosteMejorPop(j) %Si es mejor que el mejor encontrado por esa abeja
                     MejorPop(j) = beej;
                     CosteMejorPop(j) = Scorej;
                 end
            end
            
            %Actualiza valores globales
            if phij < phiMejor
                Mejor.setData(beej.data());
                phiMejor = phij;
            end
            if phij <= tolerancia && phiMejor <= tolerancia
                Scorej = cost(beej.data(),theta,n,d); %Calcula la funci�n coste de optimizaci�n principal
                 if Scorej <= CosteMejor %Si es mejor que el mejor encontrado por esa abeja
                     if   CosteMejor - Scorej >= 10e-3
                        IterSinMejora = 0;
                     end
                     Mejor.setData(beej.data());
                     CosteMejor = Scorej;
                 end
            end

            %Movimiento de las abejas
            mbee = MejorPop(j);
            Vbee = beej.dataVel(); %Velocidad de la abeja
            r = unifrnd(0,1);
            Velocidad = w*Vbee + xi*(r*(Mejor.data() - beej.data()) + (1-r)*(mbee.data() - beej.data())); %Calcula el vector de velocidad
            beej.setVel(Velocidad);
            newDato = beej.data()+Velocidad;
            for t=1:n
                if newDato(t,:) ~= zeros(1,d)
                    newDato(t,:) = newDato(t,:)/norm(newDato(t,:)); %Normaliza los vectores  
                end
            end
            beej.setData(newDato);
            Population(j) = beej; %Obtiene nueva situaci�n del individuo

        end  
        
        if w > 0
            w = w-0.1; %Se reduce el factor de inercia
        end
        
        IterSinMejora = IterSinMejora +1;
    end

end
toc

phimenor = Inf;
MejorDimension = 0;
MejorIndice=1;
for j=1:ncores
    %phiMejor{j}
    if phimenor > phiMejor{j}
        phimenor = phiMejor{j};
        MejorDimension = j+2;
        MejorIndice=j;
    end
end
ElCosteMejor = CosteMejor{MejorIndice};
ElMejor = Mejor{MejorIndice};

matlabpool close

if phimenor < 10e-6
    disp('Encontrada Representaci�n �ptima')
    dimension = MejorDimension
    best=ElMejor.data()
    ElCosteMejor
    phimenor
else
    disp('Mejor soluci�n encontrada')
    dimension = MejorDimension
    best=ElMejor.data()
    ElCosteMejor
    phimenor
end

end