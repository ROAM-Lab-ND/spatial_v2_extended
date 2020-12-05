function [qout] = quatProduct(q1, q2)
    r1 = q1(1); v1 = q1(2:4);
    r2 = q2(1); v2 = q2(2:4);
    
    rout = r1*r2 - dot(v1,v2);
    vout = r1*v2 +r2*v1 + cross(v1,v2);
    qout = [rout ; vout];
end
