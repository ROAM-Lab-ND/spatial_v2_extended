function A = randomSkew(n)
    if nargin == 0
        n = 3;
    end
    
    A = rand(n)-.5;
    A = A-A';
end

