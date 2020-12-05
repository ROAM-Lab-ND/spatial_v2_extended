function B = factorFunctions(I, v, number)
    if nargin == 2
        number = 3;
    end
    if number == 1
       B = crf(v)*I; 
    elseif number == 2
       B =  icrf( I * v) - I * crm(v); 
    else
       B = 1/2 * (crf(v)*I + icrf( I * v) - I * crm(v)); 
    end
end