    function c = cost(P,theta,n,d)
            thetaRep = CalculaTheta(n,P,d);    
            c=abs(theta-thetaRep);
    end