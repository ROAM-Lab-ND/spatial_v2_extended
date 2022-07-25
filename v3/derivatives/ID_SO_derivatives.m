function  [derivs] = ID_SO_derivatives(model, q, qd, qdd)

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
    dM_dq  =  zeros(model.NV,model.NV,model.NV);
    d2tau_dq  =  zeros(model.NV,model.NV,model.NV);
    d2tau_dqd =  zeros(model.NV,model.NV,model.NV);
    d2tau_cross =  zeros(model.NV,model.NV,model.NV);
    
for i = model.NB:-1:1
    
    for p=1:model.nv(i) % looping over each DoF of joint j

        S_p    = S{i}(:,p); Sd_p = Sd{i}(:,p); psid_p = psid{i}(:,p); psidd_p = psidd{i}(:,p);
        
        Bic_phii     = 2*factorFunctions(IC{i} ,S_p   );
        Bic_psii_dot = 2*factorFunctions(IC{i} ,psid_p);
        
        A0 = crf_bar(IC{i}*S_p);            A1 = dot(IC{i}, S_p); 
        A2 = 2*A0-Bic_phii;                 A3 = Bic_psii_dot+dot(BC{i}, S_p);
        A4 = crf_bar(BC{i}.'*S_p);          A5 = crf_bar(BC{i}*psid_p + IC{i}*psidd_p+crf(S_p)*f{i});
        A6 = crf(S_p)*IC{i}+A0;             A7 = crf_bar(BC{i}*S_p+IC{i}*(psid_p+Sd_p));
            
      ii = model.vinds{i};
      j = i;

      while j > 0
        jj = model.vinds{j}; 
        
          for t=1:model.nv(j)
              
            S_t    = S{j}(:,t);    Sd_t   = Sd{j}(:,t);
            psid_t = psid{j}(:,t); psidd_t = psidd{j}(:,t);
        
            u1 = A3.'*S_t;                                  u2 = A1.'*S_t;                               
            u3 = A3*psid_t+A1*psidd_t+A5*S_t;               u4 = A6*S_t;                                   
            u5 = A2*psid_t+A4*S_t;                          u6 = Bic_phii*psid_t + A7*S_t; 
            u7 = A3*S_t +A1*(psid_t+Sd_t);                  u8 = A4*S_t-Bic_phii.'*psid_t; 
            u9 = A0*S_t;                                    u10 = Bic_phii*S_t;  
            u11 = Bic_phii.'*S_t;                           u12 = A1*S_t;
           
               k = j;
                  while k > 0
                    kk = model.vinds{k};
                    
                     for r=1:model.nv(k)
                           S_r    = S{k}(:,r);  Sd_r   = Sd{k}(:,r); psid_r = psid{k}(:,r); psidd_r = psidd{k}(:,r);

                           p1 = u11.'*psid_r;
                           p2 = u8.'*psid_r+u9.'*psidd_r; 
                           d2tau_dq(ii(p),jj(t),kk(r)) = p2;  
                           d2tau_cross(ii(p),kk(r),jj(t)) =  -p1; 

                            if (j~=i)  

                                   d2tau_dq(jj(t),kk(r),ii(p)) =  u1.'*psid_r+ u2.'*psidd_r;     
                                   d2tau_dq(jj(t),ii(p),kk(r)) =  d2tau_dq(jj(t),kk(r),ii(p));        
                                   d2tau_cross(jj(t),kk(r),ii(p)) = p1;                               
                                   d2tau_cross(jj(t),ii(p),kk(r)) = u1.'*S_r+u2.'*(psid_r+Sd_r);                                  
                                   d2tau_dqd(jj(t),kk(r),ii(p)) =u11.'*S_r;                                         
                                   d2tau_dqd(jj(t),ii(p),kk(r)) = d2tau_dqd(jj(t),kk(r),ii(p));                
                                   dM_dq(kk(r),jj(t),ii(p)) =   S_r.'*u12;
                                   dM_dq(jj(t),kk(r),ii(p)) =  S_r.'*u12;
                            end

                             if (k~=j)
                                  d2tau_dq(ii(p),kk(r),jj(t)) = p2;  
                                  d2tau_dq(kk(r),ii(p),jj(t)) =  S_r.'*u3;                                                  
                                  d2tau_dqd(ii(p),jj(t),kk(r)) =-u11.'*S_r;     
                                  d2tau_dqd(ii(p),kk(r),jj(t)) = -u11.'*S_r;     
                                  d2tau_cross(ii(p),jj(t),kk(r)) = S_r.'*u5+u9.'*(psid_r+Sd_r);                               
                                  d2tau_cross(kk(r),jj(t),ii(p)) =  S_r.'*u6;                                 
                                  dM_dq(kk(r),ii(p),jj(t)) =  S_r.'*u9;
                                  dM_dq(ii(p),kk(r),jj(t)) =  S_r.'*u9;

                                  if(j~=i)

                                      d2tau_dq(kk(r),jj(t),ii(p)) =  d2tau_dq(kk(r),ii(p),jj(t));       
                                      d2tau_dqd(kk(r),ii(p),jj(t)) =  S_r.'*u10;      
                                      d2tau_dqd(kk(r),jj(t),ii(p)) = d2tau_dqd(kk(r),ii(p),jj(t));                 
                                      d2tau_cross(kk(r),ii(p),jj(t)) =  S_r.'*u7;                      

                                  else
                                      d2tau_dqd(kk(r),jj(t),ii(p)) =  S_r.'*u4;             
                                  end
                             else                 
                                d2tau_dqd(ii(p),jj(t),kk(r)) = -u2.'*S_r;                  

                             end

                          end
                    k  = model.parent(k);
                  end

        end 
        j  = model.parent(j);
      end
        
        
        
    end  

    if model.parent(i) > 0
        p = model.parent(i);
        IC{p} = IC{p} + IC{i};
        BC{p} = BC{p} + BC{i};
        f{p}  = f{p}  + f{i};
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

