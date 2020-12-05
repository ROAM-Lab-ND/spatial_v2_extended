function [YH, Yg, YH_rotor, Yg_rot]  = RegressorHandG(model, q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

YH = [];
Yg = [];
YH_rotor = [];
Yg_rot = [];

qd = zeros(model.NV,1);
qdd = zeros(model.NV,1);

[Yg, Yg_rot] = RegressorClassical(model,q,qd,qdd);
model.gravity = [0 0 0]';

for i = 1:model.NV
    qdd(i) = 1;
    [Y, Y_rotor] =  RegressorClassical(model,q,qd,qdd);
    qdd(i) = 0;
    YH = [ YH ; Y];
    if sum(model.NB_rot) > 0
        YH_rotor = [YH_rotor ; Y_rotor];
    end 
end
