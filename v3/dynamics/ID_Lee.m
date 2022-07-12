function  [tau, out] = ID_Lee( model, q, qd, qdd, f_tip )

% Adapted from Bryan Dongik Lee, 2018
% Original code at: https://github.com/SNURobotics/optimal-excitation/tree/master/libraries/dynamics


model = model.postProcessModel();
a_grav = model.getGravity();

if ~iscell(q)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end

v = {};
a = {};

for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  ap = model.getParentVariable(i, a, -a_grav);
  
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);

  a{i} = S{i}*qdd{i} + Xup{i}*ap + Sd{i}*qd{i};
  h{i} = model.I{i}*v{i};
end


% This line ensures that ID is symbolic or Casadi compatible
tau = q{1}(1)*0 + zeros(model.NV,1);
    
for i = model.NB:-1:1
    p = model.parent(i);  
    ii = model.vinds{i};
    if i==model.NB
        if nargin == 5
            warning('This ID does not support external forces yet');
            f_tip_transform = apply_external_forces( model.parent{end}, Xup{end}, f{end}, f_tip );
            f{i} = -f_tip_transform + model.I{i}*a{i} + crf(v{i})*h{i}; % TODO: check negative sign?
        else
            f{i} = model.I{i}*a{i} + crf(v{i})*h{i};
        end 
    else
      f{i} = Xup{i+1}.'*f{i+1} + model.I{i}*a{i} + crf(v{i})*h{i};
    end
    tau(ii) = S{i}.' * f{i};
end

out.Xup = Xup;
out.S = S;
out.Sd = Sd;
out.v = v;
out.h = h;
out.a = a;
out.f = f;
out.tau = tau;