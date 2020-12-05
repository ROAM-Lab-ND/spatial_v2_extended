function [mout] = mrpProduct(m1, m2)
    s1 = m1'*m1;
    s2 = m2'*m2;
    mout = (1-s1)*m2 + (1-s2)*m1 + 2*cross(m1,m2);
    mout = mout / (1+s1*s2 - 2*m1'*m2);
end
