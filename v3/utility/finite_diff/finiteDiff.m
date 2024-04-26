function df = finiteDiff(f, x, dx, step)
    if nargin == 3
        step = sqrt(eps);
    end
    df = (f(x+step*dx)-f(x))/step;
end