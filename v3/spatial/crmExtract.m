function  v = crmExtract( vcross )

v = zeros(6,1);
v(1:3) = skew2(vcross(1:3,1:3));
v(4:6) = skew2(vcross(4:6,1:3));
