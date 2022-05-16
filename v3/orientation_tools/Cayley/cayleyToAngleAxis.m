function angle_axis = cayleyToAngleAxis(c,n)

if nargin == 1
    n = 1;
end

angle_axis = rotToAngleAxis( cayleyToRot(c,n) );

end

