function I = symbolicInertia( body_number )
% symbolicInertia returns a symbolic Inertia matrix 
%
% Input: (optional) body_number (let's call it i)
%         
%
% Ouput: symbolic 6x6 spatial inertia

if nargin == 0
    I = inertiaVecToMat( symbolicInertialParams() );
else
	I = inertiaVecToMat( symbolicInertialParams(body_number) );
end

end

