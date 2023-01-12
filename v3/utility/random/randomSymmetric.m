function A = randomSymmetric(n)
    A = rand(n)-.5;
    A = A+A';
end

