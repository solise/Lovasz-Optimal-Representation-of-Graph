    function A = generaSolucion(n,M,d)
        %Calcula una soluci�n (pero OJO, es traspuesta)  
       	%A= rand(d,1); % Genera un vector aleatorio con d componentes.
		A= rand(d,1)+i*rand(d,1);
        for k=2:n 									% Bucle de construcci�n de la matriz A
			C=zeros(1,d);                           % Genera una matriz de ceros de dimensi�n 1xd
			for z=1:k-1 								% Bucle de construcci�n de las condiciones de la matriz
				if M(k,z)                               %Por cada k, para cada j anterior a k si son adyacentes...
					C=[C;transpose(conj(A(:,z)))];       %... A�ade a C una fila con el traspuesto-conjugado de la columna j-�sima de A
                end
            end
            nul = null(C);
			s=size(nul,2);                     % s = n�mero de columnas de la base del kernel de C (C*s=0)
			x=nul*(rand(s,1));      % Obtiene un vector del espacio lineal del n�cleo de C como combinaci�n lineal de su base 
            A=[A,x];
        end 
    end