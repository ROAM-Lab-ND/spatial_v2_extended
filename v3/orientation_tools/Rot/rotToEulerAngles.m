function [angs1, angs2] = rotToEulerAngles(R,seq)

% Gets the angle between two vectors
% assumes ax is perpendicular to a0 and a1. ax sets the sign of the result
getAng = @(a0, a1, ax) atan2( dot(ax, cross(a0,a1)) , dot(a0,a1) );

I = eye(3);
Re = @(i,th) expm(skew(I(:,i))*th);

a0 = I(:,seq(1));
a1 = a0;
b3 = R(:,seq(2));
c3 = R(:,seq(3));
c2 = c3;
    
angs2 = [0 0 0]';
for sign = [-1 1]
    angs1 = angs2;    
    b1 = sign*cross(a1,c3);
    b1 = b1/norm(b1);
    b2 = b1;
    b0 = I(:,seq(2));
    angs2(1) = getAng(b0, b1, a0);
    
    c1 = Re(seq(1), angs2(1));
    c1 = c1(:,seq(3));
    angs2(2) = getAng(c1, c2, b1);
    angs2(3) = getAng(b2, b3, c3);
end


end

