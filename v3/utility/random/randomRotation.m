function R = randomRotation()
    r = randn(3,1);
    r = r/norm(r);
    th = rand()*pi;
    R = expm(skew(r*th));
end
