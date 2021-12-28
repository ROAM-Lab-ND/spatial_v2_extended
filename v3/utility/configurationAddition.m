function new_q = configurationAddition(model,q,dq)

    USE_MCX = 0;
    
    assert(USE_MCX == 0 , 'MCX Not Yet Supported');
%     assert(1==2,'Not Yet Implemented')
    if ~isfield(model,'nq')
        model = postProcessModel(model);
    end
    if ~iscell(q)
        [q, dq] = confVecToCell(model,q,dq);
    end
    
    if USE_MCX
        new_q = MultiComplex.zeros(model.NQ,1);
    else
        new_q = q{1}(1)*0 + zeros(model.NQ,1);
    end
    
    for i = 1:model.NB
        ii = model.qinds{i};
        
        switch model.jtype{i}
            case {'Fb'}
                X = jcalc('Fb',q{i});
                v = dq{i}(4:6);
                Rup = X(1:3,1:3);
                
                % tangent element
                tang = quatR([0 ; dq{i}(1:3)/2 ]);
                    
                if ~isreal(dq{i})
                  if USE_MCX
                    qt_new = expmComplexStep(tang)* q{i}(1:4);
                  else  
                    qt_new = (eye(4) + tang)* q{i}(1:4);
                  end
                else
                  qt_new = expm(tang)* q{i}(1:4);
                end
                
                p_new  = q{i}(5:7) + Rup.'*v;
                new_q(ii) = [qt_new ; p_new];
                
            case {'S'}
                
                tang = quatR([0 ; dq{i}/2 ]);
                if ~isreal(dq{i})
                    if USE_MCX
                        new_q(ii) = expmComplexStep(tang)* q{i};
                    else
                        new_q(ii) = (eye(4) + tang)*q{i};
                    end
                else
                    new_q(ii) = expm(tang)* q{i};
                end
            case {'SO3'}
                R = reshape(q{i},[3 3]);
                so3 = -skew(dq{i});
                if ~isreal(dq{i})
                    if USE_MCX
                        new_q(ii) = reshape( expmComplexStep(so3)* R, [9 1]);
                    else
                        new_q(ii) = reshape( (eye(3) + so3)* R, [9 1]);
                    end
                else
                    new_q(ii) = reshape(expm(so3) * R, [9 1]);
                end
            case {'SE3'}
                Tup = reshape(q{i},[4 4]);
                se3 = -vecTose3(dq{i});
                if ~isreal(dq{i})
                    if USE_MCX
                        new_q(ii) = reshape( expmComplexStep(se3)*Tup, [16 1]);
                    else
                        new_q(ii) = reshape( (eye(4) + se3)*Tup, [16 1]);
                    end
                else   
                    new_q(ii) = reshape( expm(se3)*Tup, [16 1]);
                end
            otherwise
                new_q(ii) = q{i} + dq{i};
        end
    end
    if isreal(new_q)
        new_q = real(new_q);
    else if USE_MCX && new_q.order() == 1
        new_q = real(new_q) + imag(new_q)*1i;
    end
end    
