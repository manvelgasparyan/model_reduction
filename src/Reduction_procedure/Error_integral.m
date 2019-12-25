%==========================================================================
%----------------------The error integral,--------------------------------- 
%----------which gives a measure of the difference between-----------------
%---------------the original model and the reduced model-------------------
%==========================================================================
function [Error] = Error_integral(x,y) 
%It is a function of two variables x and y, which are the vectors of the
%species' concentrations in the origianl and reduced models, respectively
%--------------------------------------------------------------------------
global T; global J;
M = size(x);
n = M(1);
E_total = 0;
for j = 1:numel(J)
    p = J(j);
    for i = 1:n
        x_imp(i) = abs(y(i,p)/x(i,p) - 1); 
    end
    x_imp = x_imp';
    E(j) = Simpsons_rule(x_imp,T);
    E_total = E_total + E(j);
end

Error = E_total/(numel(J)*T);
end
%==========================================================================
