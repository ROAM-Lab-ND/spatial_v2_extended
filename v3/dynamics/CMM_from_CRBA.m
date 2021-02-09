function  [A] = CMM_from_CRBA( model, q)

% Algo from Wensing and Orin (IJHR)

assert(strcmp(model.jtype(1),'Fb'), 'First joint should be floating base');

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
qd = zeros(model.NV,1);
if ~iscell(q)
    [q,qd] = confVecToCell(model,q,qd);
end



% Compute CMM from Mass Matrix
[H] = HandC(model,q,qd);

% Joint model for floating base
[X10, Phi ] = jcalc( model.jtype{1}, q{1} );

Psi = inv(Phi);
H11 = H(1:6,1:6);

IC  = Psi'* H11 * Psi;

[~, p1G] = mcI( IC ); % extract CoM position rel to FB

R1G = X10(1:3,1:3);
X1G = [R1G  zeros(3); skew(p1G)*R1G R1G];

A = X1G' * Psi' * H(1:6, :);