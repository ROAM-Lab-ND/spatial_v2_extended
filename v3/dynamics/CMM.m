function  [A, Ad_qd] = CMM( model, q, qd)

% Algo from Orin and Goswami (AURO)
assert(~any( model.has_rotor ), 'Rotors not supported for CMM')

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

if ~iscell(q)
    if nargout == 1
        [q] = confVecToCell(model,q);
    else
        [q, qd] = confVecToCell(model, q, qd);
    end
end


A = q{1}(1)*0 + zeros(6,model.NV);
for i =1:model.NB
    IC{i} = q{1}(1)*0 + model.I{i};
end
I0 = q{1}(1)*0 + zeros(6,6);

for i = model.NB:-1:1
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i}; 
  
  if model.parent(i) ~= 0
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}.'*IC{i}*Xup{i};
  else
      I0 = I0 + Xup{i}.'*IC{i}*Xup{i};
  end
  
end

M  = I0(6,6); %Mass
pG = skew( I0(1:3,4:6)/M ); %CoM position rel. to 0
X0G = [eye(3) zeros(3) ; skew(pG) eye(3)];

for i = 1:model.NB
    p = model.parent(i);
    if p~=0
        XiG{i} = Xup{i}*XiG{p};
    else 
        XiG{i} = Xup{i}*X0G;
    end
    
    ii = model.vinds{i};
    A(:,ii) = XiG{i}.'*IC{i}*S{i};
end
%info.XiG = XiG;
%info.IC = IC;

if nargout == 2
    model_no_grav = model;
    model_no_grav.gravity = [0;0;0];
    [~, info] = ID(model_no_grav, cell2mat(q), cell2mat(qd), 0*cell2mat(qd));
    Ad_qd = X0G.'*info.f0;
end

