function q = SU2toQuat(S)
% Converts a special unitary 2x2 matrix to a quaternion

q = [ real(S(1,1) + S(2,2)) ; 
      real(S(2,1) - S(1,2)) ;
     -imag(S(2,1) + S(1,2)) ;
      imag(S(1,1) - S(2,2))]/2;
  