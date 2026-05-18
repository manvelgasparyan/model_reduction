function [time, conc] = demo_example_orig_conc (k_0, S, X_0, T, n)
    %-------------
    ds_dt = @(tau, x)original_ode_model (tau, x, k_0, S);
    %-------------
    time = 0:T/n:T;
    %-------------
    n = size(X_0, 1);
    %-------------
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.1);
    %-------------
    conc = {};
    for i = 1:n
        init_conc_i = X_0(i, :);
        [~, conc_i] = ode15s(ds_dt, time, init_conc_i, options); 
        conc{end+1} = conc_i'; 
    end
%======================
    function [ds_dt] = original_ode_model (~, s, k_0, S)
        %-----
        g = [1/(1+s(1)^2); 1/(1 + s(2)); 1];
        %-----
        v = [g(1)*k_0(1)*s(1); g(2)*k_0(2)*s(2); g(3)*k_0(3)*s(3)];
        %-----
        ds_dt = S*v; 
        %-----
    end
%======================
end
