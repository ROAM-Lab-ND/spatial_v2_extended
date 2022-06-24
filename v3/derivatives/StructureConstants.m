function [s_first_kind, s_second_kind] = StructureConstants(model, q)

    if ~isfield(model,'nq')
        model = model.postProcessModel();
    end
    if ~iscell(q)
        [q] = confVecToCell(model,q);
    end

    % S_second_kind(i,j,k) = X^i \cdot [X_j , X_k]
    % (i.e., the i-th component of the lie bracket between vector fields
    % X_j and X_k)
    s_second_kind = zeros( model.NV, model.NV, model.NV );
    
    % S_first_kind(i,j,k) = H_{i m} S_second_kind(m, j, k)
    s_first_kind  = zeros( model.NV, model.NV, model.NV );
    
    supported_joints = {'revoluteJoint','floatingBaseJoint','revoluteJointWithRotor'};
    
    
    for joint_num = 1:model.NB
       assert( any(strcmp(supported_joints, class(model.joint{joint_num}))), ...
                   ['Joint Type Unsupported In StructureConstants: ' class(model.joint{joint_num})])
        
       joint_inds =  model.vinds{joint_num} ;
       dofs = length( joint_inds );
       qj = rand( length(model.qinds{joint_num} ) ,1);

       [~, S] = model.joint{joint_num}.kinematics(model.Xtree{joint_num}, q{joint_num});
       
       [Q, ~] = qr(S);
       Psi = inv([S Q(:,dofs+1:end)])';
       Psi = Psi(:,1:dofs);

       for i = 1:dofs 
           ii = joint_inds(i);
           Si = S(:,i);
           for j = 1:dofs
               jj = joint_inds(j);
               Sj = S(:,j);
               s_second_kind(joint_inds, ii,jj ) = Psi'*crm(Si)*Sj;
           end
       end
    end
    H = HandC(model,cell2mat(q),zeros(model.NV,1));
    for i = 1:model.NV
       s_first_kind(:,:,i) = H*s_second_kind(:,:,i); % index lowering 
    end
end
