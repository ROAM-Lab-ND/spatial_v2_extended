function dot_q = configurationRates(model,q,qd)

    if ~isfield(model,'nq')
        model = postProcessModel(model);
    end
    if ~iscell(q)
        [q, qd] = confVecToCell(model,q,qd);
    end
    
    dot_q = zeros(model.NQ,1);
    
    for i = 1:model.NB
        ii = model.qinds{i};
        
        switch model.jtype{i}
            case {'Fb'}
                R = rq(q{i}(1:4));
                dot_q(ii) = [rqd(q{i}(1:4),qd{i}(1:3) ) ;
                        R * qd{i}(4:6)];
            case {'S'}
                R = rq(q{i}(1:4));
                dot_q(ii) = [rqd(q{i}(1:4),qd{i}(1:3) )];
                
            otherwise
                dot_q(ii) = qd{i};
        end
    end
end

        
