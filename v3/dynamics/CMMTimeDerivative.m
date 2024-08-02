function  [Adot, info] = CMMTimeDerivative( model, q,qd )

% Algo from Orin and Goswami (AURO)

assert(~any( model.has_rotor ), 'Rotors not supported for CMM')
%assert(strcmp(model.jtype(1),'Fb'), 'First joint should be floating base');


if ~isfield(model,'nq')
    model = postProcessModel(model);
end

if ~iscell(q)
    [q,qd] = confVecToCell(model,q,qd);
end

Adot = q{1}(1)*0 + zeros( 6, model.NV );

v = {};
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  Sd{i} = crm(v{i})*S{i};
  Id{i} = crf(v{i})*model.I{i}-model.I{i}*crm( v{i} ); % Time derivative of inertia. Derivative taken in world frame, then expressed in local.
end


A = q{1}(1)*0 + zeros(6,model.NV);
for i =1:model.NB
    IC{i} = q{1}(1)*0 + model.I{i};
    ICd{i} = q{1}(1)*0 + Id{i}; % Time derivative of composite inertia. Derivative taken in world frame, then expressed in local.
end
I0 = q{1}(1)*0 + zeros(6,6);
I0d = q{1}(1)*0 + zeros(6,6);

for i = model.NB:-1:1
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i}; 
  
  if model.parent(i) ~= 0
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}.'*IC{i}*Xup{i};
    ICd{model.parent(i)} = ICd{model.parent(i)} + Xup{i}.'*ICd{i}*Xup{i};
  else
      I0 = I0 + Xup{i}.'*IC{i}*Xup{i};
      I0d = I0d + Xup{i}.'*ICd{i}*Xup{i};
  end
end

M  = I0(6,6); %Mass
pG = skew( I0(1:3,4:6)/M ); %CoM position rel. to 0
X0G = [eye(3) zeros(3) ; skew(pG) eye(3)];
vG  = skew( I0d(1:3,4:6)/M );
vG = [zeros(3,1) ; vG];


for i = 1:model.NB
    p = model.parent(i);
    if p~=0
        XiG{i} = Xup{i}*XiG{p};
    else 
        XiG{i} = Xup{i}*X0G;
    end
    
    ii = model.vinds{i};
    A(:,ii) = XiG{i}.'*IC{i}*S{i};
    Adot(:,ii) = XiG{i}.'*(ICd{i}*S{i}+IC{i}*Sd{i});
end

% At this point, Adot has time derivative taken in world frame. Need to add
% correction for derivative taken in moving G frame.
Adot = Adot - crf(vG)*A;