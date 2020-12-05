function PrintParameterSummary_motor(model, N, M, V, C, a)
%% PrintParameterSummary_motor(model, N, M, V, C, a)
% Prints a parameter summary from the outputs of the RPNA algorithm

    for i = 1:model.NB
       fprintf('===============================\n');
        fprintf('Minimal Parameters for Body %d\n',i);
        x = 0;
        for k = 1:10
            if norm(M{i}(k,:)) > 0
                fprintf('%s\n',a{k});
                x = x+1;
            end
        end
        for k = 1:10
            if norm(M{i}(k+10,:)) > 0
                fprintf('Motor %s\n',a{k});
                x = x+1;
            end
        end
        if x == 0
            fprintf('None\n');
        end
        
        x=0;
        fprintf('\nIdentifiable Parmeters for Body %d\n',i);
        for k = 1:10
           if IsIdentifiable_rotor(model,N,i,k)
              fprintf('%s\n',a{k});
              x=x+1;
           end
        end 
        
        for k = 1:10
           if IsIdentifiable_rotor(model,N,i+model.NB,k)
              fprintf('Motor %s\n',a{k});
              x=x+1;
           end
        end 
        if x == 0
            fprintf('None\n');
        end
        

        x=0;
        fprintf('\nUnidentifiable Parmeters for Body %d\n',i);
        for k = 1:10
           if IsUnidentifiable_rotor(model,N,i,k)
              fprintf('%s\n',a{k});
              x=x+1;
           end
        end
        
        for k = 1:10
           if IsUnidentifiable_rotor(model,N,i+model.NB,k)
              fprintf('Motor %s\n',a{k});
              x=x+1;
           end
        end
        if x == 0
            fprintf('None\n');
        end
        

        fprintf('\nDim Null N(%d) = %d\n',i,20-size(N{i},1));
        fprintf('Rank VelocitySpan V(%d) = %d\n',i, rank(V{i}));  
        fprintf('Rank OuterProductSpan K(%d) = %d\n\n',i, rank(C{i}));  
        
    end