function  [dqdd_dq, dqdd_dqd,dqdd_dtau] = FD_derivatives( model, q, qd, tau )

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

qdd = FDab(model,q,qd,tau);
[q, qd, qdd] = confVecToCell(model,q,qd,qdd);


[dtau_dq, dtau_dqd] = ID_derivatives(model,q,qd,qdd);
Hinv = Hinverse(model,q);

% Relationship between FD derivatives and ID derivatives.
dqdd_dtau = Hinv;
dqdd_dq   = -Hinv*dtau_dq;
dqdd_dqd  = -Hinv*dtau_dqd;
