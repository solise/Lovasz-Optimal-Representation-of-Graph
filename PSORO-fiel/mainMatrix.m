%% PROGRAMA DE LECTURA DE FICHERO DE MATRICES DE GRAFOS PARA RO FIEL
%%  Alberto Solís Encina. Copyright 2013. 
%% NOTA: Indicar en la variable 'n' el orden de los grafos de entrada.

fid = fopen('AMqg10.txt','r');
fid2 = fopen('ResultRO10.txt', 'w');
counter=1;
n=10;
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
    sol = PSOroFiel(n, G, e);
    d = size(sol,2);
    %%Almacenar los resultados en un fichero de Texto
    
    fprintf(fid2,'%d\t%d\r\n',counter,d);
    
    counter=counter+1
    
end
fclose(fid);
fclose(fid2);