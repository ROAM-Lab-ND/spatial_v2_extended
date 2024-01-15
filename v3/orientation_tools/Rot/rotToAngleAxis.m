function aa = rotToAngleAxis(R)
aa = skew2(logm(R));
end

