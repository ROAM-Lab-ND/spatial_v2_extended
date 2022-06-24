function [s_first_kind, s_second_kind] = StructureConstants(model, q)
    
    % S_second_kind(i,j,k) = X^i \cdot [X_j , X_k]
    % (i.e., the i-th component of the lie bracket between vector fields
    % X_j and X_k)
    s_second_kind = zeros( model.NV, model.NV, model.NV );
    
    % S_first_kind(i,j,k) = H_{i m} S_second_kind(m, j, k)
    s_first_kind  = zeros( model.NV, model.NV, model.NV );
    
    for joint_num = 1:model.NB
       joint_inds =  model.vinds{joint_num} ;
       dofs = length( joint_inds );
       qj = rand( length(model.qinds{joint_num} ) ,1);

       [~, S] = jcalc(model.jtype{joint_num}, qj);
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
    H = HandC(model,q,zeros(model.NV,1));
    for i = 1:model.NV
       s_first_kind(:,:,i) = H*s_second_kind(:,:,i); % index lowering 
    end
end
