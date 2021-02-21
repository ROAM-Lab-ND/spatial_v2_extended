function a = randomConsistentInertialParams()
    J = rand(4,4)-.5;
    J = (J+J')/2;
    
    J = expm(J);
    a = pinertiaToVec(J);
end

