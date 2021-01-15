function  [o1,o2] = plux( i1, i2 )

% plux  compose/decompose Plucker coordinate transform.
% X=plux(E,r) and [E,r]=plux(X)  compose a Plucker coordinate transform X
% from its component parts E and r, and decompose it into those parts,
% respectively.  E is a 3x3 rotational coordinate transform and r is a 3D
% vector.  r is returned as a column vector, but it can be supplied as a
% row or column vector.  X is a coordinate transform corresponding to a
% shift of origin by an amount specified by r, followed by a rotation about
% the new origin as specified by E.  For example, plux(rx(1),[2 3 4]) makes
% the same transform as rotx(1)*xlt([2 3 4]).  If two arguments are
% supplied then they are assumed to be E and r, otherwise X.

if nargin == 2				% E,r --> X
    if all(size(i1) == [3 3])
        o1 = [ i1, zeros(3); -i1*skew(i2), i1 ];
    else
        o1 = [1 0 0 ; i1*[-i2(2) i2(1)]' i1];
    end
else
    X = i1;  % X --> E,r
    if all(size(X)==[6 6])			% 3D points
        E = X(1:3,1:3);
        r = -skew(E'*X(4:6,1:3));
    else					% 2D points
        E = X(2:3,2:3);
        r = [ X(2,3)*X(2,1)+X(3,3)*X(3,1); X(2,3)*X(3,1)-X(3,3)*X(2,1) ];
    end
    o1 = E;
    o2 = r;
end
