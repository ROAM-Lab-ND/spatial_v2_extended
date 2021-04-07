function  [grad_q, grad_qd, grad_tau] = modFD_derivatives( model, q, qd, tau, lambda )

% modFD derivs from modID derivs

model_no_grav = model;
model_no_grav.gravity = [0;0;0];
mu = FDab(model_no_grav,q,0*qd, lambda);

qdd = FDab(model,q,qd,tau);

[grad_q, grad_qd] = modID_derivatives(model,q,qd, qdd, -mu);
grad_tau = mu;