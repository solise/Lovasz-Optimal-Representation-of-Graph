function PSOroInt(n, M, theta, e)
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


Population = [];
MejorDimension=[]; %Guarda los mejores obtenidos en cada dimensión
for g=1:n
   MejorDimension = [MejorDimension bee(n,g)]; 
end

CosteMejor=Inf; %Coste del mejor individuo

%% Bucle principal
while (d<n) && (CosteMejor > 10e-6)

    Population = CalculaPoblacionInicEntera(n,M,d);
    MejorPop = Population; %Lista con los mejores de cada generación de cada abeja
    CosteMejorPop = ones(1,N)*Inf; %Costes de la mejor poblacion
    PhiMejorPop = ones(1,N)*Inf; %Valores de mejor distancia obtenida por individuo
    CosteMejor=Inf; %Coste del mejor individuo
    phiMejor =Inf;
    Mejor = bee(n,d);
    w=5; %Factor de inercia

    for i=1:iter %Iteraciones sobre cada dimensión
            
        for j =1:N %Para cada abeja
            beej = Population(j);
            
            phij = f(beej.data(),M,d,n); %Calcula la función de coste sobre cada individuo según las condiciones de factible

            %Actualiza Valores locales de la abeja
            if phij < PhiMejorPop(j)
                MejorPop(j) = beej;
                PhiMejorPop(j) = phij;
            end
            if phij==0 && phiMejor ==0
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
            if phij == 0 && phiMejor == 0
                Scorej = cost(beej.data(),theta,n,d); %Calcula la función coste de optimización principal
                 if Scorej <= CosteMejor %Si es mejor que el mejor encontrado por esa abeja
                    Mejor.setData(beej.data());
                     CosteMejor = Scorej;
                 end
            end

            %Movimiento de las abejas
            mbee = MejorPop(j);
            Vbee = beej.dataVel(); %Velocidad de la abeja
            Velocidad = w*Vbee + randi([1,2],1,1)*(Mejor.data() - beej.data()) + randi([1,2],1,1)*(mbee.data() - beej.data()); %Calcula el vector de velocidad
            beej.setVel(Velocidad);
            beej.setData(beej.data()+Velocidad);
            Population(j) = beej; %Obtiene nueva situación del individuo

        end  
        
        if w > 0
            w = w-1; %Se reduce el factor de inercia
        end
        
    end


MejorDimension(d) = Mejor;
d = d+1;
end


if CosteMejor < 10e-6
    disp('Encontrada Solución Optima Lovasz')
    dimension = d-1
%     EvaluaMejores = Inf*ones(1,n);
%     for g=3:n
%         MD = MejorDimension(g);
%         EvaluaMejores(g) = Evalua(MD.data(),M,g,theta,n,e);
%     end
%     [CosteMejor, IndElMejor] =  min(EvaluaMejores);
%    ElMejor = MejorDimension(IndElMejor);
    ElMejor = Mejor;
    ElMejor.data()
    CosteMejor
    phiMejor
    CosteMejorPop;
    %f(MejorPop,M,d,theta,n)
    Evalua(ElMejor.data(),M,d-1,theta,n,e)
else
    disp('Mejor solución encontrada')
    EvaluaMejores = ones(1,n)*Inf;
    for g=3:n
        MD = MejorDimension(g);
        if any(any(MD.data() ~= zeros(n,g)))
            EvaluaMejores(g)=Evalua(MD.data(),M,g,theta,n,e);
        end
    end

    [auxiliar, IndElMejor] =  min(EvaluaMejores);
    ElMejor = MejorDimension(IndElMejor);
    dimension = IndElMejor
    ElMejor.data()
    CosteMejor
    phiMejor
    %CosteMejor = cost(ElMejor.data(),theta,n,dimension)
    %phiMejor = f(ElMejor.data(),M,dimension,n)
end

%%%%%%%%%%%%%%%%%%%%%% FUNCIONES AUXILIARES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function phi = f(P,M,d,n)
          A = abs((P*P').*M);
          B = (P*P').*(eye(n)-M);  
          aux=0;
          for u=1:n
              for k=1:n
                  if u~=k
                    aux = aux + (2- rank([P(u,:);P(k,:)]));
                  end
              end
          end
          % abs(theta-thetaRep) +
          phi =  sum(sum(A)) + aux + abs( n*(n-1)/2 - e - sum(sum( (A-B) ~= 0 )) );                
    end
    
    function c = cost(P,theta,n,d)
            thetaRep = CalculaTheta(n,P,d);    
            c=abs(theta-thetaRep);
    end

    function t = CalculaTheta(n,P,d)
        %% Se calcula como el máximo de los autovalores de la suma de los proyectores
        sumProy = zeros(d,d);
        for l=1:n
            sumProy = sumProy + (P(l,:)'*P(l,:));
        end     
        t=max(eig(sumProy));
    end

     function Popul= CalculaPoblacionInicEntera(n,M,d)
        Popul=[];
        for s=1:N
            beei = bee(n,d);
            dato = randi([-1,1],n,d);
            beei.setData(dato);
            Popul = [Popul,beei];
        end
        
     end

end

