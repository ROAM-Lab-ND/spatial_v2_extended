function [FD_q, FD_qd, FD_tau, H_qq, H_qdqd, H_qdq, H_tauq] = modFD_second_derivatives( model, q, qd, tau, lambda )

% modFD derivs from modID derivs

model_no_grav = model;
model_no_grav.gravity = [0;0;0];

mu = FDab(model_no_grav,q,0*qd, -lambda);
qdd = FDab(model,q,qd,tau);

derivs = modID_second_derivatives( model, q, qd, qdd, mu);
H_qq     = derivs.dmod_dqq;
H_qdqd   = derivs.dmod_dvv;
H_qdq    = derivs.dmod_dqv';
dtau_dq  = derivs.dtau_dq;
dtau_dqd = derivs.dtau_dv;

Hinv = Hinverse(model,q);

% First Order Partials
FD_q = -Hinv*dtau_dq;
FD_qd= -Hinv*dtau_dqd;
FD_tau= Hinv;

% Second Order Partials
P = modID_derivatives_simple(model, q, mu,[FD_q.' ; FD_qd.' ; FD_tau.'].');
P_q   = P(:,1:model.NV)';
P_qd  = P(:, model.NV+1 : 2*model.NV)';
P_tau = P(:, 2*model.NV+1 : 3*model.NV)';


H_qq  = H_qq  + P_q+P_q';
H_qdq = H_qdq + P_qd;
H_tauq= P_tau;