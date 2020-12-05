function CleanA = CleanMat(A)
%% CleanA = CleanMat(A)
% set values close to 0,1, and -1 to the exact value.

inds = find(abs(A) < 1e-9); A(inds) = 0;      
inds = find(abs(A-1) < 1e-9); A(inds) = 1;
inds = find(abs(A+1) < 1e-9); A(inds) = -1;

CleanA = A;


