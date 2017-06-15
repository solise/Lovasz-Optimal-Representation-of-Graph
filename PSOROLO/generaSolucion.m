    function A = generaSolucion(n,M,d)
        %Calcula una solución (pero OJO, es traspuesta)  
       	%A= rand(d,1); % Genera un vector aleatorio con d componentes.
		A= rand(d,1)+i*rand(d,1);
        for k=2:n 									% Bucle de construcción de la matriz A
			C=zeros(1,d);                           % Genera una matriz de ceros de dimensión 1xd
			for z=1:k-1 								% Bucle de construcción de las condiciones de la matriz
				if M(k,z)                               %Por cada k, para cada j anterior a k si son adyacentes...
					C=[C;transpose(conj(A(:,z)))];       %... Añade a C una fila con el traspuesto-conjugado de la columna j-ésima de A
                end
            end
            nul = null(C);
			s=size(nul,2);                     % s = número de columnas de la base del kernel de C (C*s=0)
			x=nul*(rand(s,1));      % Obtiene un vector del espacio lineal del núcleo de C como combinación lineal de su base 
            A=[A,x];
        end 
    end