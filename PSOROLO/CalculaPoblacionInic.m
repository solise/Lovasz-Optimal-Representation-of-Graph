function Popul= CalculaPoblacionInic(N,n,M,d)
        Popul=[];
        for s=1:N
            beei = bee(n,d);
            dato = generaSolucion(n,M,d)';
            for t=1:n
              if dato(t,:) ~= zeros(1,d)
                dato(t,:) = dato(t,:)/norm(dato(t,:)); %Normaliza los vectores  
              end
            end
            beei.setData(dato);
            Popul = [Popul,beei];
        end
        
    end