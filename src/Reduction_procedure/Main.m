clear 
clc 
close
%==========================================================================
%Call the inputs of the model reduction
Model = load('HSP');
global Z B k 
%--------------------------------------------------------------------------
Z = Model.Z; B= Model.B; x_0 = Model.x_0; k = Model.k;  
Maximum_error = Model.Maximum_error;
%-------------------------------------------------------------------------- 
%There are s species, c complexes and r reactions in the network
global s c r
[s,c] = size(Z);              [~,r] = size(B); 
%==========================================================================
%Step 1: The mathematical model and the complex graph corresponding to the
%original network
%--------------------------------------------------------------------------
        %Replace al the zero elements of x_0 with a relatively small
        %number to ignore the confusion lim_{x\to 0} ln x = - \infty
        global y_0
        y_0 = x_0;
        for i =1:s
            if x_0(i) == 0
            y_0(i) = 10^(-5);
            end
        end
%--------------------------------------------------------------------------
       %The cxr outgoing matrix Delta of the original network
       global Delta
       Delta = min(B,0);
%--------------------------------------------------------------------------
       figure(1)       %The complex graphs
       subplot(1,3,1)  %The complex graph of the original network
       G = Complex_Graph(Z,B);
       plot(G, 'Layout', 'force', 'EdgeLabel', 1:r, 'EdgeColor', [0,0,1], 'LineWidth', 1, 'NodeColor', [1,0,0],  'MarkerSize', 6);
       xlabel('The complex graph corresponding to the original network' , 'Fontsize', 9, 'Fontweight', 'bold')
%==========================================================================
%Step 2: The settling time of the network
%--------------------------------------------------------------------------
        %The steady state of the network
        global T_star
        x_star = LMFsolve(@Steady_State, y_0);
        x_star = real(x_star); 
%--------------------------------------------------------------------------
        %Determination of the settling time of the network
        global opts 
        opts = odeset('AbsTol', 1e-8);
        T_star = Settling_Time(x_star);
        T_plot = max(T_star);
%==========================================================================
%Step 3: Selecting species to be eliminated from the original network
%--------------------------------------------------------------------------
        global I
        I = Species_to_Eliminate(Z);
%==========================================================================
%Step 4: The mathematical model and the complex graph corresponding to the 
         %equivalent network
%--------------------------------------------------------------------------
      global Z_eq B_eq Delta_eq s_eq c_eq r_eq
%--------------------------------------------------------------------------
       %We determine the (s_eq)x(c_eq) compex composition matrix Z_eq of
       %the equivalent network.
       Z_eq = Z;
       %Elimination of species with the set I of indices from the original 
       %network
       Z_eq(I,:) = zeros(numel(I),c); 
       %Certain complexes then become identical. So we delete all the
       %identical columns from Z_eq and keep only one of them
       Deleted_identical_complexes = [];
       Remaining_identical_complexes = [];
       for i = 1:c-1
           for j = i+1:c
               if and(Z_eq(:,i) == Z_eq(:,j), max(Z_eq(:,j))~=0)
        Deleted_identical_complexes = [Deleted_identical_complexes, j];
        Remaining_identical_complexes = [Remaining_identical_complexes, i];
               end
           end
       end
       Z_eq(:,Deleted_identical_complexes) = [];
%--------------------------------------------------------------------------
       %We determine the (c_eq)x(r_eq) incidence composition matrix B_eq of
       %the equivalent network.
       B_eq = B;
       y = numel(Remaining_identical_complexes);
       z = numel(Deleted_identical_complexes);
       for i = 1:z
           m = Deleted_identical_complexes(i);
           p = Remaining_identical_complexes(i);
           for j = 1:r
               if B_eq(p,j)==0 && B_eq(m,j)~=0
                  B_eq(p,j) = B_eq(m,j);
               end
           end
       end
       B_eq(Deleted_identical_complexes,:) = [];
%--------------------------------------------------------------------------
       %There are s_eq species, c_eq complexes and r reactions in the 
       %equivalent network
       [s_eq, c_eq] = size(Z_eq);        [~, r_eq] = size(B_eq);
%--------------------------------------------------------------------------
       %The (c_eq)x(r_eq) outgoing matrix of the equivalent network
       Delta_eq = B_eq;
       for i = 1:r_eq
           for j = 1:c_eq
               if B_eq(j,i) == 1
                  Delta_eq(j,i) = 0;
               end
           end
       end
%--------------------------------------------------------------------------
       figure(1)       %The complex graphs
       subplot(1,3,2)  %The complex graph of the equivalent network
       G_eq = Complex_Graph(Z_eq,B_eq);
       plot(G_eq, 'Layout', 'force', 'EdgeLabel', 1:r_eq, 'EdgeColor', [0,0,1], 'LineWidth', 1, 'NodeColor', [1,0,0],  'MarkerSize', 6);
       xlabel('The complex graph corresponding to the equivalent network', 'Fontsize', 9, 'Fontweight', 'bold')
%==========================================================================
%Step 5: Independent sub-networks of the equivalent network
%--------------------------------------------------------------------------        
       for i = 1:c_eq
            Index_set_eq{i} = [];
            for j = 1:s_eq
                if Z_eq(j,i) ~= 0
                   Index_set_eq{i} = [Index_set_eq{i},j];  
                end
            end
        end
%--------------------------------------------------------------------------
        %Specifies which vertice is in which connceted component 
        G_eq = Complex_Graph(Z_eq,B_eq);
        connected_comp =  conncomp(G_eq,'Type','weak');
        %Number of Linkage Classes
        [num, ~] = unique(connected_comp);
        %We define the indices of cmplexes in each linkage class
        for i = 1:numel(num)
            Linkage_class{i} = [];
            for j = 1:numel(connected_comp)
                if connected_comp(j) == num(i)
                   Linkage_class{i} = union(Linkage_class{i},j) ;
                end
            end
           % Linkage{i} = Complexes_eq(Linkage_class{i});
        end
%--------------------------------------------------------------------------
        %Determine the linkage classes that share species and combine  
        common = [];
        for i = 1:numel(num)-1
            for j = i+1:numel(num)
                L_ij = intersect(Linkage_class{i},Linkage_class{i});
                l_ij = numel(L_ij);
                if l_ij ~=0
                    Linkage_class{i} = union(Linkage_class{i},Linkage_class{j});
                    common = [common, j];
                end
            end
        end
        Linkage_class(common) = [];
%--------------------------------------------------------------------------
        %Determine the linkage classes with more than one reaction
        for i = 1:numel(Linkage_class)
            a(i) = numel(Linkage_class{i});
        end
        A = find(a>2);
        for i = 1:numel(A)
            j  = A(i);
            Linkage_more{i} = Linkage_class{j};
        end
%--------------------------------------------------------------------------
       %We define the vectors of matrices (cells) Z_link and B_link, such 
        %that their i^th elements are the complex composition and the
        %incidence matrices of the i^th linkage class, i = 1:numel(num)
        Independent_reactions = [];
        Independent_complexes = [];
        for i = 1:numel(Linkage_more)
            %Complex composition matrix
             Complex_composition_matrix{i} = Z_eq(:,Linkage_more{i});
             Z_link{i} = Complex_composition_matrix{i};
             %-------------------------------------------------------------
             %Incidence matrix
             Incidence_matrix{i} = B_eq(Linkage_more{i},:);
             B_link{i} = Incidence_matrix{i};
        end
%--------------------------------------------------------------------------
        %Determine the species in these linkage classes 
        for i = 1:numel(Linkage_more)
            Linkage_Species{i} = [];
            Z_sp = Z_link{i};
            for j = 1:s
                if nnz(Z_sp(j,:)) ~= 0
                   Linkage_Species{i} = union(Linkage_Species{i}, j);
                end
            end
        end
%==========================================================================
%Step 6: Selecting complexes to be deleted
for m = 1:numel(A)
    global Z_tilde; global B_tilde; global Delta_tilde;
%--------------------------------------------------------------------------
    %The steady state time for the Linkage class
    global T_star T
    T_linkage = T_star(Linkage_Species{i});
    T = max(T_linkage);
    t_span = [0:T/1000:T];
%--------------------------------------------------------------------------
    Z_tilde = Z_link{m};
    [~,c_tilde] = size(Z_tilde);
    B_tilde = B_link{m};
%--------------------------------------------------------------------------
        %Determine the index sets of the complexes of the linkage class
        for i = 1:c_tilde
            Index_set_tilde{i} = [];
            for j = 1:s_eq
                if Z_tilde(j,i) ~= 0
                   Index_set_tilde{i} = [Index_set_tilde{i},j];  
                end
            end
        end
%--------------------------------------------------------------------------
       %The (c_tilde)x(r_tilde) outward matrix Delta_tilde of B_tilde
       Delta_tilde = B_tilde;
       for i = 1:r_eq
           for j = 1:c_tilde
               if B_tilde(j,i) == 1
                  Delta_tilde(j,i) = 0;
               end
           end
       end
%--------------------------------------------------------------------------
        %Consider each reversible reaction as one reaction
        B_test = B_tilde;
        Reverse_reactions = [];
        for i = 1:r_eq-1
            for j = i+1:r_eq
                if B_test(:,i) == -B_test(:,j)
                   Reverse_reactions = union(Reverse_reactions,j); 
                end
            end
        end
        Rev = numel(Reverse_reactions);
        B_test(:,Reverse_reactions) = zeros(c_tilde,Rev);
%--------------------------------------------------------------------------
        %Determine the important complexes of the linkage lass
        Important_complexes = [];
        for i =1:c_tilde
            if nnz(B_test(i,:)) == 1 
               Important_complexes = [Important_complexes,i];
            end
        end
%--------------------------------------------------------------------------
        %Determine the important species of the linkage class
        global J;
        Important_species = [];
        for i = 1:numel(Important_complexes)
            p = Important_complexes(i);
           Important_species = union(Important_species,Index_set_tilde{p});
        end
%--------------------------------------------------------------------------    
        %Determine the candidate complexes for deletion from the linkage
        %class
        Original_complexes = 1:c_tilde;
        Candidates_Rao = setdiff(Original_complexes,Important_complexes);
        J = Important_species;
%-------------------------------------------------------------------------- 
        %Determine the optimal combination of complexes to be deleted from 
        %the equivalent network
        Deleted_Complexes = [];
%--------------------------------------------------------------------------
        global red;
        global del;
%--------------------------------------------------------------------------
        Threshold = Maximum_error;
        [~,x] = ode23tb(@Equivalent_Model, t_span, y_0, opts); 
        cand = Candidates_Rao;
        error = 0.1;
        red = Original_complexes;
        del = [];
%--------------------------------------------------------------------------
        while error <= Threshold
            p = numel(cand);
            red1 = red;
            del1 = del;
            error1 = error;
            for i = 1:p
                elem = [cand(i)];
                red = setdiff(red1,elem);
                del = union(del1,elem);
                [~,x1] = ode23tb(@Rao_Model, t_span, y_0, opts);
                J_error(i) = Error_integral(x,x1);
            end
            [error,z] = min(J_error);
            y = cand(z);
            red = setdiff(red1,[y]);
            del = union(del1,[y]);
            cand = setdiff(cand,[y]);
            clear J_error
        end

        Complexes_to_delete = del1;
        Complexes_to_keep = red1;
        Error(m) = error1;
        
        for i = 1:numel(Complexes_to_delete)
            j = Complexes_to_delete(i);
            for p = 1:c_eq
                if Z_eq(:,p) == Z_tilde(:,j)
                   Deleted_Complexes = union(Deleted_Complexes,p); 
                end
            end
        end
%--------------------------------------------------------------------------        
end
          Remaining_Complexes = setdiff(1:c_eq,Deleted_Complexes);
          global Deleted; global Reduced;
          Deleted = Deleted_Complexes;
          Reduced = Remaining_Complexes;
%--------------------------------------------------------------------------        
        %The complex composition matrix and the incidence matrix of the 
        %reduced model
        L = @Laplacian_Reduced;
        L_star = L(x_star);
        %Number of complexes in the reduced model
        [~,c_red] = size(L_star);
        %Number of reactions in the reduced model
        r_red = numel(find(L_star<0));
%--------------------------------------------------------------------------
        %The complex composition (s_red)x(c_red) matrix z_red of the 
        %reduced network 
        Z_red = Z_eq(:,Reduced);
%--------------------------------------------------------------------------
         %Define the vertices of the complex graph of the reduced model
         substrates_red = [];
         products_red = [];
         for i = 1:c_red
             for j = 1:c_red 
                 if L_star(i,j) < 0
                    substrates_red = [substrates_red, j]; 
                    products_red = [products_red, i];
                 end
             end
         end
%--------------------------------------------------------------------------         
         %The incidence (c_red)x(r_red) matrix B_red of the reduced network 
         B_red = zeros(c_red,r_red);
         for i = 1:r_red
             p = substrates_red(i);
             q = products_red(i);
             B_red(p,i) = -1;
             B_red(q,i) = 1;
         end
%--------------------------------------------------------------------------
         figure(1)
         subplot(1,3,3)
         G_red = Complex_Graph(Z_red,B_red);
         plot(G_red, 'Layout', 'force', 'EdgeLabel', 1:r_red, 'EdgeColor', [0,0,1], 'LineWidth', 1, 'NodeColor', [1,0,0],  'MarkerSize', 6);
         xlabel('The complex graph corresponding to the reduced network', 'Fontsize', 9, 'Fontweight', 'bold') 
%==========================================================================         
       %Plotting the concentrations
       figure(2)       
       t_span = 0:T_plot/100000:T_plot;
       [t,x] = ode23tb(@Original_Model, t_span, y_0, opts);
       [t1,y] = ode23tb(@Reduced_Model, t_span, y_0, opts);
       subplot(2,2,1)
       Species_1x = plot(t,x(:,2),  'color', [1, 0, 0], 'Linewidth', 2);
       hold on;
       Species_1y = plot(t1,y(:,2), '--', 'color', [0, 0, 1], 'Linewidth', 2);       
       hold off;
       grid on;   % Big grid
       grid minor  % Small grid
       Species_1x.Color(4) = 0.5;  % alpha
       Species_1y.Color(4) = 0.5;  % alpha
       ylabel('Concentration of X_2 (\mu mol\cdot l^{-1})' , 'Fontsize', 9, 'Fontweight', 'bold')
       legend('Original', 'Reduced')
       subplot(2,2,2)
       Species_2x = plot(t,x(:,5), 'color', [1, 0, 0], 'Linewidth', 2);
       hold on;
       Species_2y = plot(t1,y(:,5), '--', 'color', [0, 0, 1], 'Linewidth', 2);       
       hold off;
       grid on;   % Big grid
       grid minor  % Small grid
       Species_2x.Color(4) = 0.5;  % alpha
       Species_2y.Color(4) = 0.5;  % alpha       
       title({'Species'}, 'Fontsize', 9, 'Fontweight', 'bold')
       ylabel('Concentration of X_5 (\mu mol\cdot l^{-1})' , 'Fontsize', 9, 'Fontweight', 'bold')
       legend('Original', 'Reduced')
       subplot(2,2,3)
       Species_3x = plot(t,x(:,10), 'color', [1, 0, 0], 'Linewidth', 2);
       hold on;
       Species_3y = plot(t1,y(:,10), '--', 'color', [0, 0, 1], 'Linewidth', 2);       
       hold off;
       grid on;   % Big grid
       grid minor  % Small grid
       Species_3x.Color(4) = 0.5;  % alpha
       Species_3y.Color(4) = 0.5;  % alpha       
       title({'Species'}, 'Fontsize', 9, 'Fontweight', 'bold')
       ylabel('Concentration of X_{10} (\mu mol\cdot l^{-1})' , 'Fontsize', 9, 'Fontweight', 'bold')
       legend('Original', 'Reduced')
       subplot(2,2,4)
       Species_12x = plot(t,x(:,21), 'color', [1, 0, 0], 'Linewidth', 2);
       hold on;
       Species_12y = plot(t1,y(:,21), '--', 'color', [0, 0, 1], 'Linewidth', 2);       
       hold off;
       grid on;   % Big grid
       grid minor  % Small grid
       Species_12x.Color(4) = 0.5;  % alpha
       Species_12y.Color(4) = 0.5;  % alpha       
       title({'Species'}, 'Fontsize', 9, 'Fontweight', 'bold')
       ylabel('Concentration of X_{21} (\mu mol\cdot l^{-1})' , 'Fontsize', 9, 'Fontweight', 'bold')
       legend('Original', 'Reduced')         
