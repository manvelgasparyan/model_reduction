%==========================================================================
%-------The following function corresponds to the rational terms-----------
%-------in the reaction rates of the original reaction network-------------
%==========================================================================
function [d] = Denominator(x)
%--------------------------------------------------------------------------
%-------It is a function of one variables x, which is the vector of--------
%-------species' concentrations--------------------------------------------
%--------------------------------------------------------------------------
        d(1) = 1;         
        d(2) = 1;      
        d(3) = 1;     
        d(4) = 1; 
        d(5) = 1; 
        d(6) = x(7)/(1.30 + x(10) + 1.30*x(8)/0.215);  
        d(7) = x(7)/(0.215 + 0.215*x(10)/1.30 + x(8));
        d(8) = 1; 
        d(9) = 1; 
        d(10) = 1; 
        d(11) = x(7);
        d(12) = x(7); 
        d(13) = x(7); 
        d(14) = x(7);
        d(15) = 1; 
        d(16) = 1; 
        d(17) = 1; 
        d(18) = 1; 
        d(19) = 1; 
        d(20) = 1;                  
end
%==========================================================================
