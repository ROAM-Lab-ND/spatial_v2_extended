function y_V = vel_reg_Lee(V)
% assemble velocity regressor based on twist V

% Adapted from Bryan Dongik Lee, 2018
% Original code at: https://github.com/SNURobotics/optimal-excitation/tree/master/libraries/dynamics

if (length(V(:))==6 && size(V,2)==1)
    y_V = zeros(6,10);
    w = V(1:3); v = V(4:6);
    y_V(4:6,1) = v;
    y_V(:,2:4) = [-skew(v);skew(w)];
    % if Ivec = [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]:
%     y_V(1:3,5:end) = [w(1),  0  ,  0   ,  w(2),   0  , w(3); % original version from Lee's thesis
%                        0  , w(2),  0   ,  w(1),  w(3),   0 ;
%                        0  ,  0  , w(3) ,   0  ,  w(2), w(1) ];
    % if Ivec = [Ixx, Iyy, Izz, Izy, Ixz, Iyx]:
    y_V(1:3,5:end) = [w(1),  0  ,  0   ,   0  ,  w(3), w(2); % modified to match Patrick's ordering of parameters
                       0  , w(2),  0   ,  w(3),   0  , w(1);
                       0  ,  0  , w(3) ,  w(2),  w(1),  0   ];

elseif mod(length(V(:)),6)==0 && size(V,2)==1
    i = 1;
    j = 1;
    while i < length(V(:))
        ii = i:i+5; % i.e. 1:6
        jj = j:j+9; % i.e. 1:10
        y_V(ii,jj) = vel_reg_Lee(V(ii));
        i = i+6;
        j = j+10;
    end

else
    warning("Wrong size for velocity regressor.")
    y_V = zeros(6,10);

end

end