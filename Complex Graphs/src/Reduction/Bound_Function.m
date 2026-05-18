%==========================================================================
%----Incidence function which is defined equal to one in the 3% bound of---
%-----------------x_star, and 0 out of that bound--------------------------
%==========================================================================
function [incidence_function] = Bound_Function(x, y)
         if x>=0.97*y && x<=1.03*y
            incidence_function = 1;
         else
            incidence_function = 0;
         end
end
%==========================================================================
