function mat = icrf(f)

ii = 1;
while ii < length(f)
    inds = ii:ii+5;
    ff = f(inds);
    mat(inds,inds) = [-skew(ff(1:3)) -skew(ff(4:6));-skew(ff(4:6)) zeros(3)];
    ii = ii+6;
end
