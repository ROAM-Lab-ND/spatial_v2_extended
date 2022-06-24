function [Gamma, out] = Christoffel(model,q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q] = confVecToCell(model,q);
end
supported_joints = {'revoluteJoint','floatingBaseJoint','revoluteJointWithRotor'};


for i = 1:model.NB
    [Xup{i}, S{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i});
    IC{i} = q{1}(1)*0 + model.I{i}; 
    S{i}  = S{i};
    IC{i} = IC{i};
end

for k = model.NB:-1:1
    assert( any(strcmp(supported_joints, class(model.joint{k}))), ['Joint Type Unsupported In Christoffel: ' class(model.joint{k})])  
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
               
               f1 = Xup{i}.'* f1;
               f2 = Xup{i}.'* f2;

               i = model.parent(i);
           end
           B = Xup{j}.'*B*Xup{j};
           j = model.parent(j);
        end
    end
    
    if model.parent(k) > 0
        p = model.parent(k);
        IC{p}= IC{p} + Xup{k}.'* IC{k} * Xup{k};
    end
end

out.S  = S;
out.IC = IC;