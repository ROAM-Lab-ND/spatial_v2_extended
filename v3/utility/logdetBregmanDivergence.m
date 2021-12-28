function d = logdetBregmanDivergence(P,Q)
    d = -log( det(P) / det(Q)) + trace(Q\P) - length(P);
end