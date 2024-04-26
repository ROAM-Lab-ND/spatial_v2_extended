function B = rightJacobianOfExp(ad_phi)
    B = eye(size(ad_phi,1));
    for i = 1:40
        B = B + ad_phi^i *(-1)^i / factorial(i+1);
    end
end

