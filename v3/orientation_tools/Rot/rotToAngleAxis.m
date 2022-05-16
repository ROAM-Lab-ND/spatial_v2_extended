function aa = rotToAngleAxis(R)
aa = skew(logm(R));
end

