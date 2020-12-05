

c = rand(3,1);
R = cayleyToRot(c);
wR = rand(3,1);
wL = R*wR;
mrp = rotToMRP(R);
a = rotToAngleAxis(R);
q = rotToQuat(R);


dt = sqrt(eps);
Rnew = R*angleAxisToRot(wR*dt);

%% Check Quat Rates
qnew = rotToQuat(Rnew);
qdot_numer = (qnew-q)/dt;
qdot  = quatRateRight(q,wR);
qdot2 = quatRateLeft(q,wL);

assert( isEqual(qdot, qdot_numer, 1e-6 ) )
assert( isEqual(qdot2, qdot_numer, 1e-6 ) )

%% Check Cayley Rates
cnew = rotToCayley(Rnew);
cdot_numer = (cnew-c)/dt;

cdot  = cayleyRateRight(c,wR);
cdot2 = cayleyRateLeft(c,wL);

assert( isEqual(cdot, cdot_numer, 1e-6 ) );
assert( isEqual(cdot2, cdot_numer, 1e-6 ) );

%% Check MRP Rates
mnew = rotToMRP(Rnew);
mdot_numer = (mnew - mrp)/dt;

mdot  = mrpRateRight(mrp, wR);
mdot2 = mrpRateLeft(mrp, wL);

assert( isEqual(mdot, mdot_numer, 1e-6 ) )
assert( isEqual(mdot2, mdot_numer, 1e-6 ) )

%% Check Angle Axis Rates
anew = rotToAngleAxis(Rnew);
adot_numer = (anew-a)/dt;

adot = angleAxisRateRight(a, wR);
adot2 = angleAxisRateLeft(a, wL);

assert( isEqual(adot, adot_numer, 1e-6 ) )
assert( isEqual(adot2, adot_numer, 1e-6 ) )







