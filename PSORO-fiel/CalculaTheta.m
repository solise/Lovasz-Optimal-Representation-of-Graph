
    function t = CalculaTheta(n,P,d)
        %% Se calcula como el m�ximo de los autovalores de la suma de los proyectores
        sumProy = zeros(d,d);
        for l=1:n
            sumProy = sumProy + (P(l,:)'*P(l,:));
        end     
        t=max(eig(sumProy));
    end