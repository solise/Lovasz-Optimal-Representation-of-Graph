function phi = f(P,M,d,n,e)
          A = abs((P*P').*M);
          B = (P*P').*(ones(n,n)-M);  
          aux=0;
          for u=1:n
              for k=1:n
                  if u~=k
                    aux = aux + (2- rank([P(u,:);P(k,:)]));
                  end
              end
          end
          % abs(theta-thetaRep) +
%           uno = sum(sum(A))
%           dos = aux
%           tres = abs( n*(n-1)/2 - e - sum(sum( (A-B) ~= 0 )) )
%             (abs(A-B) > 10e-16)-eye(n)
%           laA = A
%           laB = B
%           (abs((A-B)) > 10e-16) - eye(n)
          phi =  sum(sum(A)) + aux + abs( n*(n-1)/2 - e - sum(sum( (abs(A-B) > 10e-8)-eye(n) ))/2 );                
    end