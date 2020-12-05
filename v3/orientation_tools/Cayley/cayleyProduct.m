function [cout] = cayleyProduct(c1, c2)
    cout = (c1+c2+cross(c1,c2))/(1-c1'*c2);
end
