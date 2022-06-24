function [YH, Yg]  = RegressorHandG(model, q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

YH = [];

qd = zeros(model.NV,1);
qdd = zeros(model.NV,1);

Yg = RegressorClassical(model,q,qd,qdd);
model.gravity = [0 0 0]';

for i = 1:model.NV
    qdd(i) = 1;
    Y =  RegressorClassical(model,q,qd,qdd);
    qdd(i) = 0;
    YH = [ YH ; Y];
end
