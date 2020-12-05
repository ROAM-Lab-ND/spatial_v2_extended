function PrintParameterSummary(model, N, M, V, C, a)
%% PrintParameterSummary(model, N, M, V, C, a)
% Prints a parameter summary from the outputs of the RPNA algorithm

    for i = 1:model.NB
        
        fprintf('===============================\n');
        fprintf('Minimal Parameters for Body %d\n',i);
        x=0;
        for k = 1:10
            if norm(M{i}(k,:)) > 0
                fprintf('%s\n',a{k});
                x=x+1;
            end
        end
        if x == 0
            fprintf('None\n');
        end

        fprintf('\nIdentifiable Parmeters for Body %d\n',i);
        x=0;
        for k = 1:10
           if IsIdentifiable(model,N,i,k)
              fprintf('%s\n',a{k});
              x=x+1;
           end
        end
        if x == 0
            fprintf('None\n');
        end


        fprintf('\nUnidentifiable Parmeters for Body %d\n',i);
        x=0;
        for k = 1:10
           if IsUnidentifiable(model,N,i,k)
              fprintf('%s\n',a{k});
              x=x+1;
           end
        end
        if x == 0
            fprintf('None\n');
        end

        fprintf('\nDim Null N(%d) = %d\n',i,10-size(N{i},1));
        fprintf('Dim VelocitySpan V(%d) = %d\n',i, rank(V{i}));  
        fprintf('Dim OuterProductSpan K(%d) = %d\n\n',i, rank(C{i}));  
        
    end