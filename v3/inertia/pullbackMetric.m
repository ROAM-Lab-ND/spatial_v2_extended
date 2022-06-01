function [M] = pullbackMetric(params)
    M = zeros(10,10);
    P = inertiaVecToPinertia(params);
    P_inv = inv(P);
    v_i = zeros(10,1);
    v_j = zeros(10,1);
    for i = 1 : 10
       for j = 1 : 10
           v_i = zeros(10,1); v_j = zeros(10,1);
           v_i(i) =1 ; v_j(j) = 1;
           V_i = inertiaVecToPinertia(v_i); V_j = inertiaVecToPinertia(v_j);
           M(i,j) = trace(P_inv*V_i*P_inv*V_j);
       end
    end
end