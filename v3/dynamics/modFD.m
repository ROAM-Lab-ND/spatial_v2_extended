function  out = modFD( model, q, qd, tau, mu )

qdd = FDab(model,q,qd,tau);
out = mu.'*qdd;
