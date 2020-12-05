function q = mrpToQuat(m)
mm = m'*m;
q = [ (1-mm) ; 2*m ]/(1+mm);
end

