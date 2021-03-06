%--------------------------------------------------------------------------
%----The following function evaluates the integral of a given function f--- 
%-------------over the interval [0,b] using the Simpson's rule------------- 
function I = Simpsons_rule(f,b)
    n = numel(f)-1;
    h = b/n; 
    I = h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end-1))+f(end));
end
%--------------------------------------------------------------------------
%-----------------------------END------------------------------------------
%--------------------------------------------------------------------------
