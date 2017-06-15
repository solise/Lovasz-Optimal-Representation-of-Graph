function d=CS_graph(G)
%%
%% generate the orthogonality representation (columns of resulting matrix) of a graph given by G
%%
% ************************


e=0.00000000000001;								% Precisicón numérica
tol = 0.01;                                     % Tolerancia del producto escalar
n=size(G,1);                                    % número de vértices = nº de filas de la matriz (cuadrada)
d=3; 											% Asume que existe un nodo conectado con otros 2
bad=true;
iter=1;
index=0;

%t1=clock;

while bad   %  Hacer mientras bad -> mientras no sea buena la representación
		A=rand(d,1)+i*rand(d,1); % Genera un vector aleatorio complejo con d componentes.
		bad=false;
		for k=2:n 									% Bucle de construcción de la matriz A
			C=zeros(1,d);                           % Genera una matriz de ceros de dimensión 1xd
			for j=1:k-1 								% Bucle de construcción de las condiciones de la matriz
				if G(k,j)                               %Por cada k, para cada j anterior a k si son adyacentes...
					C=[C;transpose(conj(A(:,j)))];       %... Añade a C una fila con el traspuesto-conjugado de la columna j-ésima de A
                end
            end
            nul = null(C);
			s=size(nul,2);                     % s = número de columnas de la base del kernel de C (C*s=0)
			x=nul*(rand(s,1)+i*rand(s,1));      % Obtiene un vector del espacio lineal del núcleo de C como combinación lineal de su base
			for j=1:k-1 								% Bucle de comprobación
				if ((rank([A(:,j),x],tol)<2) || ((abs(A(:,j)'*x)<e) && G(j,k)==0)) 	% comprueba paralelismo u ortogonalidad adicional
					bad=true;
					break;
                end
			end
			if (bad && (index>=iter))
				d=d+1;
                if d>n
					bad=false;
					d=n;
                end
                index=0;
				break
            elseif (bad && (index<iter))
                index=index+1;
                break
			else
				A=[A,x]; 							% Comprobación pasada, se añade x a la siguiente columna de A
            end
		end
end
A;

	dimension=d;					%print dimension
	%B = conj(transpose(A))*A;				% Genera la matriz de Gramm B desde la matriz A
	%for k=1:size(B,1)					    % Cambia la matriz de Gramm B a una matriz ortogonal (bucle que recorre toda la matriz B)
	%	for j=1:size(B,2)
	%		B(k,j)=abs(B(k,j));
	%		if abs(B(k,j))<e
	%			B(k,j)=1;
	%		else
	%			B(k,j)=0;
    %        end
	%	end
    %end
    
    % La matriz B resultante debe corresponder a la matriz de adyacencia
    % para que el resultado sea correcto (en caso contrario se tienen
    % adyacencias adicionales).

	%additional_orthogonalities=ones(1,n)*(B-G)*ones(n,1)/2	%count ones in difference B-G, divide by 2 (additional orthogonalities)
	%B-G % El conjunto de vectores es correcto si dicha matriz resulta nula
    %tiempo = etime(clock,t1)
end


