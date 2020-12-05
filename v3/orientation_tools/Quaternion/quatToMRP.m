function c = quatToMRP(q)
if q(1) < -.99
    q = -q;
end
c = q(2:4)/(1+q(1));
end

