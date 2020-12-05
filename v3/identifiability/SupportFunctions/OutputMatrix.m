function C = OutputMatrix(L1, L2)
%% C = OutputMatrix(L1, L2)
% Gives the output matrix C(L1, L2) described in the accompanying paper
% each column of L1 (L2) is a vector l1 (l2)
% each row of C, cj, satisfies cj'*a = l1'*I*l2 for some l1, l2
% where a and I are the paried parameter vector and matrix form of any spatial inertia

    if all( size(L1) == size(L2) )
        if norm(L1-L2,'fro') < eps
            L1 = RangeBasis(L1);
            L2 = L1;
        else
            L1 = RangeBasis(L1);
            L2 = RangeBasis(L2);
        end
    else
        L1 = RangeBasis(L1);
        L2 = RangeBasis(L2);
    end
    C = OutputMatrix_NoReduce(L1,L2);
    C = RangeBasis(C')';
end

function C = OutputMatrix_NoReduce(L1, L2)
    C = zeros(0,10);
    for i = 1:size(L1,2)
        li = L1(:,i);
        for j = 1:size(L2,2)
            lj = L2(:,j);
            C = [C; Output_InnerProduct(li,lj)];
        end
    end
end

