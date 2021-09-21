function [f, g, H] = entropicDivergence(p, q)
    P = inertiaVecToPinertia(p);
    Q = inertiaVecToPinertia(q);
    f = logdetBregmanDivergence( P , Q );
    if nargout > 1
        iQ = inv(Q);
        iP = inv(P);
        G = iQ-iP;
        H = zeros(10,10);
        g = zeros(10,1);
        for i = 1:10
            ei = zeros(10,1);
            ei(i) = 1;
            Ji = inertiaVecToPinertia(ei);
            g(i) = trace( Ji * G );
            if nargout > 2
                for j=i:10
                    ej = zeros(10,1);
                    ej(j) = 1;
                    Jj = inertiaVecToPinertia(ej);
                    H(i,j) = trace(iP*Ji*iP*Jj);
                    H(j,i) = H(i,j);
                end
            end
        end
    end
end