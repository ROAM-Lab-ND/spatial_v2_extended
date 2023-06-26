function  [derivs] = ID_SO_derivatives(model, q, qd, qdd )

% SO Partial Derivatives of Inverse Dynamics for multi-DoF joints
% Contributors - Shubham Singh, singh281@utexas.edu
%              - Dr. Patrick Wensing, pwensing@nd.edu 

crf_bar = @(x)[-skew(x(1:3)), -skew(x(4:6)) ; -skew(x(4:6)) ,zeros(3,3)]; % cross product matrix for force vectors with a bar over top

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end
if sum(model.has_rotor) > 1
    error('ID_derivatives does not support rotors');
end


a_grav = get_gravity(model);
IC = model.I;
I = model.I;

% Forward sweep
for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
    Xup{i} = XJ * model.Xtree{i};
    if model.parent(i) == 0
        v{i}  = zeros(6,1);
        a{i}  = -a_grav;
        Xup0{i} = Xup{i};
    else
        Xup0{i} = Xup{i}*Xup0{model.parent(i)};
        v{i}  = v{model.parent(i)};
        a{i}  = a{model.parent(i)};
    end
    Xdown0{i} = inv(Xup0{i});
    S{i} = Xdown0{i}*S{i};
    vJ{i}= S{i}*qd{i};
    aJ{i} = crm(v{i})*vJ{i} + S{i}*qdd{i};
    psid{i} = crm(v{i})*S{i};
    psidd{i}= crm(a{i})*S{i} + crm(v{i})*psid{i};    
    v{i} = v{i} + vJ{i};
    a{i} = a{i} + aJ{i};
    IC{i} = Xup0{i}.'*I{i}*Xup0{i};
    Sd{i} = crm(v{i})*S{i};
    BC{i} = 2*factorFunctions(IC{i},v{i});
    f{i}  =  IC{i}*a{i} + crf(v{i})*IC{i}*v{i};
end

% Backward sweep: Can be parallelized across branches
for j = model.NB:-1:1
    if model.parent(j) > 0
        p = model.parent(j);
        IC{p} = IC{p} + IC{j};
        BC{p} = BC{p} + BC{j};
        f{p}  = f{p}  + f{j};
    end
end

% Can be parallelized across all j,d
for j = model.NB:-1:1
    for d=1:model.nv(j) % looping over each DoF of joint j
        dd = model.vinds{j}(d);
        S_d    = S{j}(:,d);
        Sd_d   = Sd{j}(:,d);
        psid_d = psid{j}(:,d);
        psidd_d = psidd{j}(:,d);
        Bic_phii     = 2*factorFunctions(IC{j} ,S_d   );
        Bic_psii_dot = 2*factorFunctions(IC{j} ,psid_d);
        
        A1 = dot(IC{j}, S_d); 
        A2 = Bic_psii_dot+dot(BC{j}, S_d); 
        A3 = crf_bar(IC{j}*S_d);
        
        % These temps are 6 x n (where n =#dofs)
        T1(:,dd) = IC{j}*S_d;  
        T2(:,dd) = - BC{j}.'*S_d;
        T3(:,dd) = BC{j}*psid_d + IC{j}*psidd_d+crf_bar(f{j})*S_d;
        T4(:,dd) = BC{j}*S_d+IC{j}*(psid_d+Sd_d);
        
        % These temps are 36 x n (where n =#dofs)
        D1(:,dd) = A1(:);
        D2(:,dd) = A2(:);
        D3(:,dd)  = Bic_phii(:);
        D4(:,dd)  = A3(:);
    end 
end

dM_dq  =  zeros(model.NV,model.NV,model.NV);
d2tau_dq  =  zeros(model.NV,model.NV,model.NV);
d2tau_dqd =  zeros(model.NV,model.NV,model.NV);
d2tau_cross =  zeros(model.NV,model.NV,model.NV);
   
% Can be parallelized over all j,d,k,c
for j = model.NB:-1:1
    jj = model.vinds{j};
    st_j = model.subtree_vinds{j};
    succ_j = model.successor_vinds{j};

    for d=1:model.nv(j) % looping over each DoF of joint j
        k = j;
        dd = model.vinds{j}(d);
        S_d = S{j}(:,d);
        Sd_d = Sd{j}(:,d);
        psid_d = psid{j}(:,d);
        psidd_d = psidd{j}(:,d);
        
        for k = model.ancestors{j} % looping over all ancestors of joint j (including joint j)
            
            for c=1:model.nv(k)    % looping over each DoF of joint k
                cc = model.vinds{k}(c);
                S_c    = S{k}(:,c);
                Sd_c   = Sd{k}(:,c);
                psid_c = psid{k}(:,c);
                
                % these temps are 36 x 1
                t1 = S_d    * psid_c.'    ;    t1 = t1(:);        
                t2 = S_d    * S_c.'       ;    t2 = t2(:);       
                t3 = psid_d * psid_c.'    ;    t3 = t3(:);        
                t4 = S_d    * psidd{k}(:,c).'; t4 = t4(:);        
                t5 = S_d    * (Sd_c+psid_c).'; t5 = t5(:);
                t8 = S_c*S_d.'    ; t8 = t8(:);  
                
                % these temps are 6 x 1
                p1 = crm( psid_c)*S_d;  
                p2 = crm(psidd{k}(:,c)) * S_d;
                
                d2tau_dq(st_j,dd,cc) = -t3.'*D3(:,st_j) -p1.'*T2(:,st_j) + p2.'*T1(:,st_j);     
                d2tau_cross(st_j,dd,cc) =-t1.'*D3(:,st_j);                      
                
                if (k<j)
                    t6 = S_c*psid_d.' ; t6 = t6(:); 
                    t7 = S_c*psidd_d.'; t7 = t7(:);  
                 
                    p3 = crm(S_c)*S_d;   
                    p4 = crm(Sd_c+psid_c)*S_d-2*crm(psid_d)*S_c;                    
                    p5 = crm(S_d)*S_c; 
                                          
                    d2tau_dq (st_j,cc,dd) =  d2tau_dq(st_j,dd,cc);  
                    d2tau_dqd(st_j,cc,dd) = -t2.'*D3(:,st_j);    
                    d2tau_dqd (st_j,dd,cc) =  d2tau_dqd(st_j,cc,dd);                                 
                    d2tau_cross(st_j,cc,dd) = -t6.'*D3(:,st_j)-p3.'*T2(:,st_j)+p4.'*T1(:,st_j);    
                    d2tau_dq  (cc,st_j,dd) = t6.'*D2(:,st_j) + t7.'*D1(:,st_j) - p5.'*T3(:,st_j);                                                                                           
                    d2tau_cross(cc,st_j,dd) = t6.'*D3(:,st_j) - p5.'*T4(:,st_j);                    
                    d2tau_dqd(cc,jj,dd) = (S_d'*IC{j}*crm(S_c)+S_c'*crf(S_d)*IC{j})*S{j}; 
                    dM_dq(cc,st_j,dd) = t8.'*D4(:,st_j);
                    dM_dq(st_j,cc,dd) = dM_dq(cc,st_j,dd);
                    
                    if (isempty(succ_j))==0
                        t9 = S_c*(Sd_d+psid_d).'; t9=t9(:);
                        d2tau_dqd (cc,succ_j,dd) =  t8.'*D3(:,succ_j);      
                        d2tau_cross(cc,dd,succ_j) =  t8.'*D2(:,succ_j) +t9.'*D1(:,succ_j) ;                      
                        d2tau_dq(cc,dd,succ_j) =  d2tau_dq(cc,succ_j,dd);
                        d2tau_dqd(cc,dd,succ_j) = d2tau_dqd(cc,succ_j,dd);              
                    end
                end
                if (isempty(succ_j))==0
                    d2tau_dq  (dd,cc,succ_j) = t1.'*D2(:,succ_j) + t4.'*D1(:,succ_j);    
                    d2tau_dqd (dd,cc,succ_j) = t2.'*D3(:,succ_j);                                 
                    d2tau_cross(dd,succ_j,cc) = t1.'*D3(:,succ_j);                                         
                    d2tau_dq  (dd,succ_j,cc) =  d2tau_dq (dd,cc,succ_j);        
                    d2tau_dqd (dd,succ_j,cc) =  d2tau_dqd(dd,cc,succ_j);                                  
                    d2tau_cross(dd,cc,succ_j) =  t2.'*D2(:,succ_j)+t5.'*D1(:,succ_j); 
                    dM_dq(cc,dd,succ_j) = t8.'*D1(:,succ_j);
                    dM_dq(dd,cc,succ_j) = dM_dq(cc,dd,succ_j);
                end
                if(k==j)
                    d2tau_dqd(st_j,dd,cc) = -t2.'*D1(:,st_j);             
                end
            end
            k  = model.parent(k);
        end
    end

end

derivs.d2tau_dq=d2tau_dq;
derivs.d2tau_dqd=d2tau_dqd;
derivs.d2tau_cross=d2tau_cross;
derivs.dM_dq=dM_dq;

end

function Idot = dot(I, v)
    Idot = crf(v)*I - I*crm(v);
end

