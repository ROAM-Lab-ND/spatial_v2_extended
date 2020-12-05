function [q_cell, v_cell, vd_cell, vd2_cell] = confVecToCell(model,q,v,vd, vd2)
    
    qj = 0;
    vj = 0;
    
    if ~isfield(model,'nq')
        model.nq = ones(model.NB,1);
        model.nv = ones(model.NB,1);
    end
    
    q_cell = cell(model.NB,1);
    v_cell = cell(model.NB,1);
    vd_cell = cell(model.NB,1);
    vd2_cell = cell(model.NB,1);
    
    for i = 1:model.NB
        
        q_cell{i} = q(qj+1 : qj+model.nq(i) );
        if nargin > 2
            v_cell{i} = v( vj+1 : vj+model.nv(i) );
        end
        if nargin > 3
            vd_cell{i} = vd( vj+1 : vj+model.nv(i) );
        end
        if nargin > 4
            vd2_cell{i} = vd2( vj+1 : vj+model.nv(i) );
        end
        
        qj = qj+model.nq(i);
        vj = vj+model.nv(i);
    end
    
end