function d = affineInvariantDistance( P, Q )
    d = norm( 1/sqrt(2) * log( eig(P\Q) ) );
end

