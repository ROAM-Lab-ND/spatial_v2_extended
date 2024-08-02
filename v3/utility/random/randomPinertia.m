function J = randomPinertia()
    J = rand(4,4)-.5;
    J = (J+J');
    J = expm(J);
end

