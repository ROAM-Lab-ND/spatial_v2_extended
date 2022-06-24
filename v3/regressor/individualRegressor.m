function Y = individualRegressor(a,v,w,factorFunction)

if nargin <4
    factorFunction = @ (I,v) factorFunctions(I,v);
end

if nargin < 3
    w = v;
end

Y = zeros(6,10);
for i =1:10
   ei = zeros(10,1);
   ei(i) = 1;
   I = inertiaVecToMat(ei);
   Y(:,i) = I*a + factorFunction(I,v)*w;
end
