%% PROGRAMA DE LECTURA DE FICHERO DE GRAFOS PARA USAR ALGORITMO DE GRAMM
%%  Alberto Solís Encina. Copyright 2012. Máster en matemática computacional
%% NOTA: Indicar en la variable 'n' el orden de los grafos de entrada.

fid = fopen('ForbidenGraphs5.txt','r');
fid2 = fopen('ResultFG5.txt', 'w');
counter=1;
n=5;
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
    d = CS_graph(G);
    %%Almacenar los resultados en un fichero de Texto
    
    fprintf(fid2,'%d\t%d\r\n',counter,d);
    
    counter=counter+1
    
end
fclose(fid);
fclose(fid2);