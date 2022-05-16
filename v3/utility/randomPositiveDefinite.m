function A = randomPositiveDefinite(n)
    A = expm( randomSymmetric(n) );
end

