function new_q = configurationAddition(model,q,dq)

%     assert(1==2,'Not Yet Implemented')
    if ~isfield(model,'nq')
        model = postProcessModel(model);
    end
    if ~iscell(q)
        [q, dq] = confVecToCell(model,q,dq);
    end
    
    new_q = zeros(model.NQ,1);
    
    for i = 1:model.NB
        ii = model.qinds{i};
        
        switch model.jtype{i}
            case {'Fb'}
                X = jcalc('Fb',q{i});
                v = dq{i}(4:6);
                Rup = X(1:3,1:3);
                
                dquat = angleAxisToQuat( dq{i}(1:3) );
                qt_new = quatProduct(q{i}(1:4), dquat );
                p_new  = q{i}(5:7) + Rup'*v;
                
                new_q(ii) = [qt_new ; p_new];
                
            case {'S'}
                dquat = angleAxisToQuat( dq{i} );
                new_q(ii) = quatProduct(q{i}(1:4), dquat );
            case {'SO3'}
                R = reshape(q{i},[3 3]);
                so3 = -skew(dq{i});
                if norm( imag(dq{i})) > 0
                    new_q(ii) = reshape( so3* R + R, [9 1]);
                else
                    new_q(ii) = reshape(expm(so3) * R, [9 1]);
                end
            case {'SE3'}
                Tup = reshape(q{i},[4 4]);
                se3 = vecTose3(dq{i});
                if norm( imag(dq{i})) > 0
                    % Special case for complex step on SE3
                    new_q(ii) = reshape( -se3*Tup + Tup, [16 1]);
                else   
                    new_q(ii) = reshape( expm(-se3)*Tup, [16 1]);
                end
            otherwise
                new_q(ii) = q{i} + dq{i};
        end
    end
end

        
