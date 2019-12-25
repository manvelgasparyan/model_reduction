%==========================================================================
%--------The following function is used to determine the settling----------
%--------time of a given biochemical reaction network----------------------
%==========================================================================
function [T_star] = Settling_Time(x_star)
%--------------------------------------------------------------------------
%-------It is a function of one variables x_star, which is the steady------
%-------state of a given biochemcial network network-----------------------
%--------------------------------------------------------------------------
        global s y_0 opts
        T_star = zeros(s,1); %the output vector
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
        non_vanishing_species = [];
        X_star = round(x_star);
        for i = 1:s
            if X_star(i) ~=0
               non_vanishing_species = [non_vanishing_species, i]; 
            end
        end
        nvs = numel(non_vanishing_species);
        vanishing_species = setdiff(1:s,non_vanishing_species);
        vs = numel(vanishing_species);
%-------------------------------------------------------------------------- 
        %The steady state time of the original network
        t_test = Steady_State(y_0);
        t_test = abs(t_test);
        t_test = max(t_test);
        x_test = abs(x_star - y_0);
        x_test = max(x_test);
        dt = x_test/t_test;
        T1 = 10^3*dt;
        t_span = 0:dt:T1;
        global opts;
        opts = odeset('AbsTol', 1e-8);
        [~,x] = ode23tb(@Original_Model, t_span, y_0, opts);
        for i = 1:nvs
            j = non_vanishing_species(i);
            while Bound_Function(x(end,j), x_star(j)) == 0
                T1 =  T1 + T1;
                break
            end
        end
        t_span = 0:dt:T1;
        [t,x] = ode23tb(@Original_Model, t_span, y_0, opts);
%-------------------------------------------------------------------------- 
       %non vanishing species
       t_star = zeros(s,1);
        for i = 1:nvs
            j = non_vanishing_species(i);
            position = [];
            y = x(:,j);
            z = x_star(j);   
            for p = 1:numel(y)-1
                if and(Bound_Function(y(p),z) == 0,Bound_Function(y(p+1),z) == 1)
                   position = [position p+1];
                end
            end
            t_1 = max(position);
            t_star(i) = t(t_1);
        end
        T_star = t_star;
end
%==========================================================================
