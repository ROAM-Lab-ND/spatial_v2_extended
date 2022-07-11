function Y = individualRegressor(a,v,w,factorFunction)

if nargin <4
    factorFunction = @ (I,v) factorFunctions(I,v);
end

if nargin < 3
    w = v;
end

% Hotfix for symbolic REMOVE
Y = zeros(6,10);
%symY = cell(1,10); % REMOVE
for i =1:10
   ei = zeros(10,1);
   ei(i) = 1;
   I = inertiaVecToMat(ei);
   Y(:,i) = I*a + factorFunction(I,v)*w;
   %symY{i} = I*a + factorFunction(I,v)*w; %REMOVE
end
%Y = [symY{1} symY{2} symY{3} symY{4} symY{5} symY{6} symY{7} symY{8} symY{9} symY{10}]; %REMOVE
