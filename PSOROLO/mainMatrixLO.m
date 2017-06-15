%% PROGRAMA DE LECTURA DE FICHERO DE MATRICES DE GRAFOS PARA RO FIEL
%%  Alberto Solís Encina. Copyright 2013. 
%% NOTA: Indicar en la variable 'n' el orden de los grafos de entrada.

fid = fopen('AMqg9.txt','r');
fid2 = fopen('ResultROLO9.txt', 'w');
FicheroDatos = importdata('Quantum9.txt');
Thetas = FicheroDatos.data(:,2);
counter=1;
n=9;

while ~feof(fid)
    G=zeros(n,n);
    for i=1:n
        fscanf(fid,'%c',1);
        A = fscanf(fid, '%d', n);
        if ~(isempty(A))
            G(i,:)=A;
        end
    end
    if any(any(G))
        G;
    else 
        break;
    end
    fscanf(fid,'%c',2);
    fscanf(fid,'\n',1);
    
    %%Realizar operaciones sobre el grafo G
    e=sum(sum(G))/2;
    tic
    [sol, cost, d] = PSOro1Dimension(n, G, Thetas(counter), e);
    t=toc
    %%Almacenar los resultados en un fichero de Texto
    
    fprintf(fid2,'%d\t%d\t%f\t%f\r\n',counter,d, cost,t);
    
    counter=counter+1
    
end
fclose(fid);
fclose(fid2);