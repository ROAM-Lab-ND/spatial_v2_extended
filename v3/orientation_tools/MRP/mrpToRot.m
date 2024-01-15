function R = mrpToRot(m)
mm = m'*m;
R = eye(3) + (8*skew2(m)^2+4*(1-mm)*skew2(m))/(1+mm)^2;
end

