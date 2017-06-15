   function costes = Evalua(P,M,d,theta,n,e)
          costes = Inf;
          A = abs ((P*P').*M);
          B = (P*P').*(ones(n,n)-M);  
         thetaRep = CalculaTheta(n,P,d);
          aux=0;
          for u=1:n
              for k=1:n
                  if u~=k
                    aux = aux + (2- rank([P(u,:);P(k,:)]));
                  end
              end
          end
          % abs(theta-thetaRep) +
          costes =  abs(theta-thetaRep) + sum(sum(A)) + aux + abs( n*(n-1)/2 - e - sum(sum( (abs(A-B) > 10e-16)-eye(n) ))/2 );
          
          
          function t = CalculaTheta(n,P,d)
        %% Se calcula como el máximo de los autovalores de la suma de los proyectores
        sumProy = zeros(d,d);
        for l=1:n
            sumProy = sumProy + (P(l,:)'*P(l,:));
        end     
        t=max(eig(sumProy));
          end 
          
    end