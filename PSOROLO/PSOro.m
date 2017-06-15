function [Sol,CosteMejor,dimension]=PSOro(n, M, theta, e)
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

iter = 1000; %Número de iteraciones
N = 50; %Tamaño del enjambre
d=3; %Dimension
tolerancia=10e-15;


Population = [];
MejorDimension=[]; %Guarda los mejores obtenidos en cada dimensión
for g=1:n
   MejorDimension = [MejorDimension bee(n,g)]; 
end

CosteMejor=Inf; %Coste del mejor individuo

%% Bucle principal
while (d<n) && (CosteMejor > 10e-3)
  
    Population = CalculaPoblacionInic(N,n,M,d);
    MejorPop = Population; %Lista con los mejores de cada generación de cada abeja
    CosteMejorPop = ones(1,N)*Inf; %Costes de la mejor poblacion
    PhiMejorPop = ones(1,N)*Inf; %Valores de mejor distancia obtenida por individuo
    CosteMejor=Inf; %Coste del mejor individuo
    phiMejor =Inf;
    Mejor = bee(n,d);
    w=0.6; %Factor de inercia
    xi = 1.7; %Factor de contracción
    IterSinMejora = 0;

    while CosteMejor > 10e-2 && IterSinMejora <= iter        
        for j =1:N %Para cada abeja
            beej = Population(j);
            
            phij = f(beej.data(),M,d,n,e); %Calcula la función de coste sobre cada individuo según las condiciones de factible

            %Actualiza Valores locales de la abeja
            if phij < PhiMejorPop(j)
                MejorPop(j) = beej;
                PhiMejorPop(j) = phij;
            end
            if phij<=tolerancia && phiMejor <= tolerancia
                 Scorej = cost(beej.data(),theta,n,d); %Calcula la función coste de optimización principal
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
                Scorej = cost(beej.data(),theta,n,d); %Calcula la función coste de optimización principal
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
            beej.setData(newDato/norm(newDato));
            Population(j) = beej; %Obtiene nueva situación del individuo

        end  
        
        if w > 0
            w = w-0.1; %Se reduce el factor de inercia
        end
        
        IterSinMejora = IterSinMejora +1;
    end

Mejor
MejorDimension(d) = Mejor;
d = d+1;
end


if CosteMejor <= 10e-2
    disp('Encontrada Solución Optima Lovasz')
    %dimension = d-1
    dimension = d
    ElMejor = Mejor;
    Sol=ElMejor.data()
    CosteMejor
    phiMejor
    %f(MejorPop,M,d,theta,n)
    %Evalua(ElMejor.data(),M,d-1,theta,n,e);
else
    disp('No encontrada representacion LO')
    EvaluaMejores = ones(1,n)*Inf;
    for g=3:n
        MD = MejorDimension(g);
        %if MD ~= zeros(n,g)
            EvaluaMejores(g)=Evalua(MD.data(),M,g,theta,n,e);
        %end
    end
    [CosteMejor, IndElMejor] =  min(EvaluaMejores);
    ElMejor = MejorDimension(IndElMejor);
    if CosteMejor <= 10e-1
        %dimension = IndElMejor
        dimension=d
        Sol=ElMejor.data()
    else
        dimension = 0
        Sol=[]
    end
    CosteMejor
    phiMejor
end


end