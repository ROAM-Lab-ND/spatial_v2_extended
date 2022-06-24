function [Gamma, out] = Christoffel(model,q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q] = confVecToCell(model,q);
end
if sum(model.has_rotor) > 1
    error('Christoffel does not support rotors (yet)');
end

for i = 1:model.NB
    IC{i} = q{1}(1)*0 + model.I{i};
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
    Xup{i} = XJ * model.Xtree{i};
    if model.parent(i) == 0
        Xup0{i} = Xup{i};
    else
        Xup0{i} = Xup{i}*Xup0{ model.parent(i) };
    end
    
    S{i} = Xup0{i}\S{i};
    IC{i} =  Xup0{i}.'*IC{i}*Xup0{i};
end

for k = model.NB:-1:1
    for k_ind = 1:length( model.vinds{k} )
        kk = model.vinds{k}(k_ind);
        sk = S{k}(:,k_ind);
        B = factorFunctions(IC{k},sk);
    
        j = k;
        while j > 0
           jj = model.vinds{j};
           f1 = B  * S{j};
           f2 = B.'* S{j};
           i=j;
           while i > 0
               ii = model.vinds{i};
               
               Gamma(ii,jj,kk) = S{i}.'*f1;  
               if j < k
                Gamma(ii,kk,jj) = S{i}.'*f1;
               end

               Gamma(jj,ii,kk) = ( S{i}.'*f2 ).'; 
               if i < k
                Gamma(jj,kk,ii) = ( S{i}.'*f2 ).'; 
               end

               Gamma(kk,ii,jj) =  -S{i}.'*f2;
               if i < j
                Gamma(kk,jj,ii) = ( -S{i}.'*f2).';
               end
               
               i = model.parent(i);
           end
           j = model.parent(j);
        end
    end
    
    if model.parent(k) > 0
        p = model.parent(k);
        IC{p}= IC{p} + IC{k};
    end
end

out.S  = S;
out.IC = IC;