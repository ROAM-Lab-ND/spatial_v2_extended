clear
cayley = rand(3,1);
quat = cayleyToQuat(cayley);
R = cayleyToRot(cayley);
angle_axis = rotToAngleAxis(R);
rpy = quatToRpy(quat);
mrp = quatToMRP(quat);

[rpy, rpy2] = rotToFixedAngles(R,[1 2 3]); % Two Euler Angle Solutions

%% Check Functions that convert to quaternions
q2 = rotToQuat(R);
q3 = cayleyToQuat(cayley);
q4 = angleAxisToQuat(angle_axis);
q5 = rpyToQuat(rpy);
q6 = mrpToQuat(mrp);

assert( isEqualModSign(quat, q2) )
assert( isEqualModSign(quat, q3) )
assert( isEqualModSign(quat, q4) )
assert( isEqualModSign(quat, q5) )
assert( isEqualModSign(quat, q6) )


%% Check Functions that convert to rotations
R2 = quatToRot(quat);
R3 = cayleyToRot(cayley);
R4 = angleAxisToRot(angle_axis);
R5 = rpyToRot(rpy);
R6 = mrpToRot(mrp);

assert( isEqual(R, R2) )
assert( isEqual(R, R3) )
assert( isEqual(R, R4) )
assert( isEqual(R, R5) )
assert( isEqual(R, R6) )


%% Check Functions that convert to cayley
c2 = quatToCayley(quat);
c3 = rotToCayley(R);
c4 = angleAxisToCayley(angle_axis);
c5 = rpyToCayley(rpy);
c6 = mrpToCayley(mrp);

assert( isEqual(cayley, c2) )
assert( isEqual(cayley, c3) )
assert( isEqual(cayley, c4) )
assert( isEqual(cayley, c5) )
assert( isEqual(cayley, c6) )


%% Check Functions that convert to angle axis
a2 = quatToAngleAxis(quat);
a3 = rotToAngleAxis(R);
a4 = cayleyToAngleAxis(cayley);
a5 = rpyToAngleAxis(rpy);
a6 = mrpToAngleAxis(mrp);


assert( isEqual(angle_axis, a2) )
assert( isEqual(angle_axis, a3) )
assert( isEqual(angle_axis, a4) )
assert( isEqual(angle_axis, a5) )
assert( isEqual(angle_axis, a6) )


%% Check Functions that convert to rpy
r2 = quatToRpy(quat);
r3 = rotToRpy(R);
r4 = cayleyToRpy(cayley);
r5 = angleAxisToRpy(angle_axis);
r6 = mrpToRpy(mrp);


assert( isEqual(rpy, r2) | isEqual(rpy2, r2) )
assert( isEqual(rpy, r3) | isEqual(rpy2, r3) )
assert( isEqual(rpy, r4) | isEqual(rpy2, r4) )
assert( isEqual(rpy, r5) | isEqual(rpy2, r5) )
assert( isEqual(rpy, r6) | isEqual(rpy2, r6) )


%% Check funcitons that convert to mrp
m2 = quatToMRP(quat);
m3 = rotToMRP(R);
m4 = cayleyToMRP(cayley);
m5 = angleAxisToMRP(angle_axis);
m6 = rpyToMRP(rpy);


assert( isEqual(mrp, m2) )
assert( isEqual(mrp, m3) )
assert( isEqual(mrp, m4) )
assert( isEqual(mrp, m5) )
assert( isEqual(mrp, m6) )


%% Now Check quaternion functions

quat1 = quat;
quat2 = rand(4,1);
quat2 = quat2/norm(quat2);

result = quatProduct(quat1, quat2);

result2 = quatL(quat1)*quat2;
result3 = quatR(quat2)*quat1;

assert(isEqual(result, result2))
assert(isEqual(result, result3))

R1 = quatToRot(quat1);
R2 = quatToRot(quat2);

Rout = R1*R2;
R3   = quatToRot(result);
assert(isEqual(Rout, R3))

%% Now Check Cayley functions
c1 = cayley;
c2 = rand(3,1);
cout = cayleyProduct(c1,c2);

R1 = cayleyToRot(c1);
R2 = cayleyToRot(c2);
Rout = R1*R2;

c3 = rotToCayley(Rout);
R3 = cayleyToRot(cout);

assert(isEqual(cout, c3))
assert(isEqual(Rout, R3))


%% Now Check mrp functions

shadow = mrpShadow(mrp);
Rshadow = mrpToRot(shadow);

assert( isEqual(Rshadow, R) )

m1 =mrp;
m2 = rand(3,1);

mout = mrpProduct(m1,m2);

R1 = mrpToRot(m1);
R2 = mrpToRot(m2);
R3 = R1*R2;
m3 = rotToMRP(R3);

Rout = mrpToRot(mout);

assert( isEqual(Rout, R3) );
assert( isEqual(mout, m3) | isEqual( mout, mrpShadow(m3)) );
