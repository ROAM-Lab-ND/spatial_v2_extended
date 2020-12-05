function R = mrpToRot(m)
mm = m'*m;
R = eye(3) + (8*skew(m)^2+4*(1-mm)*skew(m))/(1+mm)^2;
end

