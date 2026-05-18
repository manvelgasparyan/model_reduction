import os, libsbml, sys, numpy, time, networkx, sympy, json, re
from sympy import simplify, symbols, sympify, nsimplify
from tkinter import filedialog, Tk
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from tkinter import filedialog, Tk
import matplotlib.pyplot as plt
from sympy import symbols
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sympy import *
#=======================================================================================  
#======================================================================================= 
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"           # Reaction ID
#-----------------------------------------------------------------
numpy.set_printoptions(threshold=numpy.inf) # Display the full array
#-----------------------------------------------------------------
os.system('cls') # Clear the terminal before output
#=======================================================================================  
#=======================================================================================
#=====================MAIN FUNCTIONS====================================================
#======================================================================================= 
"""
This function uses a graphical file dialog to allow the user to select one or more files, checks if any files are selected, 
and then returns both the full paths and the base file names of the selected files. If no files are selected, the program exits.
"""
#----------------------------------------------------------------------------
def get_file():
    #-----------
    root = Tk()
    root.withdraw()  
    #-----------
    file_path = filedialog.askopenfilename(
        parent=root,
        title="Select SBML file",
        filetypes=[("All files", "*.*")]
    )
    #-----------
    root.destroy()
    #-----------
    if not file_path:
        print(f"\n{BOLD}No file selected.\n")
        sys.exit(1)
    #-----------
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    #-----------
    return file_path, file_name
#=================================================================================================================================
#=================================================================================================================================
"""
This function loads an SBML file using a string input, file_path, which specifies the file's location. It returns a tuple containing 
an instance of libsbml.Model, representing the extracted biological model, and an instance of libsbml.SBMLDocument, which represents 
the entire SBML document, including any metadata.
"""
#----------------------------------------------------------------------------
def get_model(file_path: str) -> tuple[libsbml.Model, libsbml.SBMLDocument]:
    #---
    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromFile(file_path)
    #---
    errors = document.getNumErrors()
    if errors > 0:
        print(f"{BOLD}Error loading the SBML file:{RESET}")
        for e in range(errors):
            print(f"\t{document.getError(e).getMessage()}")
        #sys.exit(1)
    #---
    model = document.getModel()
    #---
    if model is None:
        print(f"{BOLD}No model found in the SBML file.\n{RESET}")
        sys.exit(1)
    #---
    return model, document
#=================================================================================================================================
#=================================================================================================================================
"""
This function takes a model instance as input. It returns a tuple containing four elements: a list of symbolic representations of 
non-constant species, their initial concentrations, a list of symbolic representations of constant species, and their corresponding 
initial values.
"""
#----------------------------------------------------------------------------
def model_species (model):
    #----------------------------
    all_species = [species.getId() for species in model.getListOfSpecies()]
    #----------------------------
    constant_species = list({species.getId() for species in model.getListOfSpecies() if species.getConstant() or species.getBoundaryCondition()})
    constant_species_sym = symbols(constant_species)
    constant_species_indices = [all_species.index(string) for string in constant_species]
    #----------------------------
    species = [item for item in all_species if item not in constant_species]
    species_sym = symbols(species)
    species_indices = [all_species.index(string) for string in species]
    #----------------------------
    full_initial_values = [string.getInitialAmount() if string.isSetInitialAmount() else string.getInitialConcentration()
                            for string in (model.getSpecies(i) for i in range(model.getNumSpecies()))]
    #----------------------------
    constant_species_values = [full_initial_values[i] for i in constant_species_indices]
    initial_conc = [full_initial_values[i] for i in species_indices]
    #----------------------------
    return species_sym, initial_conc, constant_species_sym, constant_species_values
#=================================================================================================================================
#=================================================================================================================================
"""
This function takes a model as input and constructs a stoichiometric matrix that represents the relationships between metabolites 
and reactions. The output is a NumPy array where rows correspond to metabolites and columns correspond to reactions
"""
#----------------------------------------------------------------------------
def stoichiometric_matrix (model):
    #---
    reactions = model.getListOfReactions()
    metabolites = model.getListOfSpecies()
    metabolite_indices = {met.getId(): i for i, met in enumerate(metabolites)}
    #----
    #Initiate the stoich matrix
    stoich = numpy.zeros((len(metabolites), len(reactions)))
    #----
    for i, reaction in enumerate(reactions):
        for reactant in reaction.getListOfReactants():
            metabolite_id = reactant.getSpecies()
            coeff = reactant.getStoichiometry()
            stoich[metabolite_indices[metabolite_id], i] -= coeff
        for product in reaction.getListOfProducts():
            metabolite_id = product.getSpecies()
            coeff = product.getStoichiometry()
            stoich[metabolite_indices[metabolite_id], i] += coeff
    #--- 
    buffered_species_id = [species.getId() for species in model.getListOfSpecies() if species.getConstant()]
    boundary_species_id = [species.getId() for species in model.getListOfSpecies() if species.getBoundaryCondition()]
    #---
    constant_species_id = list(dict.fromkeys(boundary_species_id + buffered_species_id))
    #---
    species_id = [species.getId() for species in model.getListOfSpecies()]
    buffered_indices_list = [species_id.index(string) for string in constant_species_id]
    #---
    stoich = numpy.delete(stoich, buffered_indices_list, axis=0)
    #---
    reversible_reactions = [j for j in range(model.getNumReactions()) if model.getReaction(j).getReversible()]
    #---
    new_columns = {j + 1: numpy.where(-stoich[:, j] == -0, 0, -stoich[:, j]) for j in reversible_reactions }    
    #---
    insert_positions = sorted(new_columns.keys(), reverse=True)
    #---
    stoich_new = stoich.copy()
    for pos in insert_positions:
        stoich_new = numpy.insert(stoich_new, pos, new_columns[pos], axis=1)
    #---
    return stoich_new
#=================================================================================================================================
#=================================================================================================================================
"""
This function takes a model as input and retrieves the global parameters and kinetic parameters associated with the reactions in 
the model. It outputs two lists: the first list contains the IDs of all parameters, and the second list contains their values.
"""
#----------------------------------------------------------------------------
def model_parameters(model):
    #----------------------------
    # Get global parameters and their values
    parameters = model.getListOfParameters()
    parameter_id = [param.getId() for param in parameters]
    parameter_values = [param.getValue() for param in parameters]
    #----------------------------
    # Get kinetic parameters and their values from reaction kinetics
    kinetic_parameters = [
        (param.getId(), param.getValue())
        for reaction in model.getListOfReactions()
        if reaction.getKineticLaw() is not None
        for param in reaction.getKineticLaw().getListOfParameters()
    ]
    #----------------------------
    # Unzip kinetic parameters into separate lists
    kinetic_parameter_id, kinetic_parameter_values = zip(*kinetic_parameters) if kinetic_parameters else ([], [])
    #----------------------------
    # Combine global and kinetic parameters
    parameter_id.extend(kinetic_parameter_id)
    parameter_values.extend(kinetic_parameter_values)
    #----------------------------
    parameter_sym = symbols(parameter_id)
    #----------------------------
    return parameter_sym, parameter_values
#=================================================================================================================================
#=================================================================================================================================
"""
This function takes a model as input and retrieves the IDs and sizes of its compartments. It returns a tuple containing symbolic 
representations of the compartment IDs and a list of their corresponding sizes.
"""
#----------------------------------------------------------------------------
def model_compartments(model):
    #----------------------------
    # Get compartment IDs and sizes
    compartments = model.getListOfCompartments()
    compartment_id = [compartment.getId() for compartment in compartments]
    compartment_size = [compartment.getSize() for compartment in compartments]
    #----------------------------
    compartment_sym = symbols(compartment_id)
    #----------------------------
    return compartment_sym, compartment_size
#=================================================================================================================================
#=================================================================================================================================
"""
This function takes a model as input, extracting kinetic rate laws for its reactions. It retrieves function definitions, arguments, 
and their expressions, and then substitutes the functions in the reaction rate formulas with their respective expressions. 
The output is a symbolic list of modified reaction rate formulas with substituted function definitions.
"""
#----------------------------------------------------------------------------
def model_rates (model):
    #==========================
    def model_functions (model):
        #---
        def extract_arguments(Inside):
            Arguments = []

            #Function to check if a character is a mathematical symbol
            def is_math_symbol(char):
                return char in ('+', '-', '*', '/', '^', '(', ')')

            #Loop until a mathematical symbol is encountered
            while Inside:
                #Find the index of the first comma
                comma_index = Inside.find(',')

                #If there is no comma or the comma is at the start, break the loop
                if comma_index == -1 or comma_index == 0:
                    break

                #Extract the part before the first comma
                part = Inside[:comma_index].strip()

                #If the part contains a mathematical symbol, break the loop
                if any(is_math_symbol(c) for c in part):
                    break

                #Add the part to Arguments and remove it from Inside along with the comma
                Arguments.append(part)
                Inside = Inside[comma_index+1:].strip()
            return Arguments, Inside
        #---
        #Get the list of functions in the model
        function_id_list = []
        function_list = []
        #---
        for i in range(model.getNumFunctionDefinitions()):
            function = model.getFunctionDefinition(i)
            function_id = function.getId()
            math_ast = function.getMath()
            math_string = libsbml.formulaToString(math_ast)
            function_id_list.append(function_id)
            function_list.append(math_string)
        #---
        arguments = []
        expressions = []
        for i in range(len(function_id_list)):
            #---
            inside_parentheses_i = re.search(r'lambda\((.*)\)', function_list[i]).group(1)
            #print("Inside parantheses:", inside_parentheses_i)
            #---
            arguments_i, expression_i = extract_arguments(inside_parentheses_i)
            #---
            arguments.append(arguments_i)
            expressions.append(expression_i)
            #---
        return function_id_list, function_list, arguments, expressions
    #==============================
    reaction_rate_species_ID = []
    for reaction in model.getListOfReactions():
        reaction_rate_formula = reaction.getKineticLaw().getFormula() if reaction.getKineticLaw() else None
        if reaction_rate_formula:
           reaction_rate_species_ID.append(reaction_rate_formula)  
    #---
    function_id, _, arguments, expressions = model_functions(model)
    #---
    reaction_rate_species_ID_new = reaction_rate_species_ID.copy()
    #---
    for j in range(len(reaction_rate_species_ID_new)):
        #---
        f_j =  [[func, re.search(re.escape(func) + r'\((.*?)\)', reaction_rate_species_ID_new[j]).group(1)] for func in function_id if re.search(re.escape(func) + r'\((.*?)\)', reaction_rate_species_ID_new[j])]
        for i in range(len(f_j)):
            elem_ji = f_j[i]
            #---
            new_id = elem_ji[0]
            new_args = elem_ji[1].split(',')
            new_name = new_id+ '(' + ','.join(new_args) + ')'
            #---
            old_arg = arguments[function_id.index(new_id)]
            old_expres = expressions[function_id.index(new_id)]
            #---
            new_expess = old_expres
            for z in range(len(old_arg)):
                new_expess = new_expess.replace(old_arg[z], new_args[z]) 
            #---
            reaction_rate_species_ID_new[j] = reaction_rate_species_ID_new[j].replace(new_name, new_expess)     
    #---
    compartment_sym, _ = model_compartments(model)
    compartment_list = [str(vector) for vector in compartment_sym]
    #---
    parameter_sym, _ = model_parameters(model)
    parameter_list = [str(vector) for vector in parameter_sym]
    #---
    species_sym, _, constant_species_sym, _ = model_species (model)
    species_list = [str(vector) for vector in species_sym]
    constant_species_list = [str(vector) for vector in constant_species_sym]
    #---
    variables = {name: symbol for name, symbol in zip(species_list + constant_species_list + parameter_list + compartment_list, species_sym + constant_species_sym + parameter_sym + compartment_sym)}
    rates_sym_rev  = [sympy.expand(item) for item in sympify(reaction_rate_species_ID_new, locals=variables)]
    #---
    reversible_reactions = [j for j in range(model.getNumReactions()) if model.getReaction(j).getReversible()]
    #---
    rates_sym_new = rates_sym_rev.copy()
    #---
    for j in reversible_reactions:
        rates_sym_new_j = rates_sym_rev[j]
        forward_rate = sympy.Add(*[term for term in sympy.Add.make_args(rates_sym_new_j) if term.as_coeff_Mul()[0] > 0])
        rates_sym_new[j] = forward_rate
    #---
    new_columns = {j + 1: - sympy.Add(*[term for term in sympy.Add.make_args(rates_sym_rev[j]) if term.as_coeff_Mul()[0] < 0]) for j in reversible_reactions }    
    insert_positions = sorted(new_columns.keys(), reverse=True)
    #---
    for pos in insert_positions:
        rates_sym_new = numpy.insert(rates_sym_new, pos, new_columns[pos] )
    #---
    rates_sym_final = [item for item in rates_sym_new]
    #---
    return rates_sym_final
#=================================================================================================================================
#=================================================================================================================================
"""
This function takes a model as input, extracting kinetic rate laws for its reactions. It retrieves function definitions, arguments, 
and their expressions, and then substitutes the functions in the reaction rate formulas with their respective expressions. 
The output is a symbolic list of modified reaction rate formulas with substituted function definitions.
"""
#----------------------------------------------------------------------------
def stoichiometric_ode_model (stoichiometric, constant_species_sym, constant_species_values, rates_sym, parameter_sym, parameter_values, compartment_sym):
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    S = sympy.Matrix(stoichiometric)
    #---
    v1 = sympy.Matrix([expr.subs({ci: 1 for ci in compartment_sym}) for expr in rates_sym.copy()])
    #---
    ds_dt_sym = list(S*v1)
    #---
    v2 =  sympy.Matrix([expr.subs({parameter_sym[i]: parameter_values[i] for i in range(len(parameter_sym))}) for expr in v1.copy()])
    #---
    v =  sympy.Matrix([expr.subs({constant_species_sym[i]: constant_species_values[i] for i in range(len(constant_species_sym))}) for expr in v2.copy()])
    #---
    ds_dt = S*v
    #------------------------------------------------------------------
    #------------------------------------------------------------------  
    return ds_dt_sym, ds_dt
#=================================================================================================================================
#=================================================================================================================================
"""
This function computes the steady-state concentrations of species in a biochemical reaction system, along with their corresponding 
times to reach steady state. 
"""
#----------------------------------------------------------------------------
def equilibrium (ds_dt, T, n, species_sym, initial_conc, delta, method_int):
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    def fun_int(t, y):
        #---
        subs_dict = dict(zip(species_sym, y))
        #---
        dydt = [float(sympy.re(expr.subs(subs_dict))) for expr in ds_dt]
        #---
        return dydt
    #---------------
    def bound_function(x, y):
        if (1-delta)* y <= x <= (1+delta) * y:
            return 1
        else:
            return 0
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    solutions = solve_ivp(fun_int, (0, T), initial_conc, t_eval=numpy.linspace(0, T, n), method=method_int)#, atol=1e-5, rtol=1e-5) 
    solutions_time = solutions.t
    solutions_concentrations = solutions.y
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    steady_state = solutions_concentrations[:, -1] 
    #------------------------------------------------------------------            
    steady_state_list = steady_state.tolist()    
    #------------------------------------------------------------------
    tau = []
    for i in range(len(species_sym)):
        bound_i = [bound_function(solutions_concentrations[i,j],steady_state[i]) for j in range(len(solutions_time))]
        index = []
        for j1 in range(len(solutions_time)-1):
                if bound_i[j1] == 0 and bound_i[j1+1] == 1:
                    index.append(j1)
        id = max(index)
        tau_i = solutions_time[id]
        tau.append(tau_i)            
    #------------------------------------------------------------------ 
    return steady_state_list, tau
#=================================================================================================================================
#=================================================================================================================================
"""
This function calculates the concentrations of species over a specified interval by integrating the system of ordinary differential equations 
"""
#----------------------------------------------------------------------------
def ode_model_integrations (ds_dt, T, n, species_sym, species_initial_values, method_int, file_name, verbose):
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    def fun_int(t, y):
        #---
        subs_dict = dict(zip(species_sym, y))
        #---
        dydt = [float(sympy.re(expr.subs(subs_dict))) for expr in ds_dt]
        #---
        return dydt
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    solutions_final = solve_ivp(fun_int,  (0, T), species_initial_values, t_eval=numpy.linspace(0, T, n), method=method_int)
    solution_concentrations = solutions_final.y
    solution_time = solutions_final.t
    #------------------------------------------------------------------
    if verbose:
        #------------
        fig = plt.figure()
        fig.canvas.manager.set_window_title(f"Concentration plots: {file_name}")
        #------------
        # Plot the solution
        for i in range(len(species_sym)):
            plt.plot(solution_time, solution_concentrations[i], label=f'{species_sym[i]}')
        plt.xlabel('Time')
        plt.ylabel('Concentration')
        plt.title('Solution of the System of Differential Equations')
        plt.legend()
        #------------
        plt.show(block=False) 
    #------------------------------------------------------------------         
    return solution_time, solution_concentrations 
#=================================================================================================================================
#=================================================================================================================================
"""
This function computes the symmetrized error integral between concentrations in two different model.
"""
#----------------------------------------------------------------------------
def simmetrized_error_integral (x, y, tau):
    #------------------------
    def simpsons_rule (f, b):
        #---
        n = len(f) - 1
        #---
        h = b/n
        #--- 
        I = (h/3)*(f[0] + 2*numpy.sum(f[2:n-2:2]) + 4*numpy.sum(f[1:n-1:2]) + f[n-1])
        #---
        return I
    #------------------------
    def error_integral (f_1, f_2, b):
        #---
        theta = 0
        #---
        for i in range(len(f_1)):
            f_i = [numpy.abs(f_2[i,j]/f_1[i,j] - 1) for j in range(len(f_1[0])) if f_1[i,j] != 0]
            b_i = b[i]
            if f_i:
               theta += simpsons_rule(f_i,b_i)/b_i
        #---
        theta_final = theta/len(f_1)  
        #---
        return theta_final
    #------------------------
    theta_1 = error_integral (x, y, tau)
    #---
    theta_2 = error_integral (y, x, tau)
    #---
    Theta = (theta_1 + theta_2)/2
    #---
    return Theta
#=================================================================================================================================
#=================================================================================================================================
"""
This function computes the species-reaction graph of the model. 
"""
#----------------------------------------------------------------------------
def species_reaction_graph(stoichiometric, species_sym, file_name, verbose=0):
    #------------
    graph = networkx.DiGraph()
    #------------
    #Label the species and the reaction nodes
    species = [str(vector) for vector in species_sym]
    reactions = [f"R{i}" for i in range(1, len(stoichiometric[0]) + 1)]
    #-----------
    graph.add_nodes_from(species)  
    graph.add_nodes_from(reactions) 
    #------------
    edges = []
    #------------
    for i in range(len(stoichiometric)):
        for j in range(len(stoichiometric[0])):
            if stoichiometric[i][j] < 0:
                edges.append((species[i], reactions[j]))
    #------------
    for i in range(len(stoichiometric)):
        for j in range(len(stoichiometric[0])):
            if stoichiometric[i][j] > 0:
               edges.append((reactions[j], species[i]))
    #------------
    graph.add_edges_from(edges)
    #------------
    #====================
    if verbose:
       #------------
       fig = plt.figure()
       fig.canvas.manager.set_window_title(f"Species-reaction graph: {file_name}")
       #------------
       pos = networkx.spring_layout(graph)
       #pos = networkx.circular_layout(graph)
       #------------
       networkx.draw_networkx_nodes(graph, pos, nodelist=species, node_color='lightblue', node_size=1000, alpha=1.0, 
                                   node_shape='o',  edgecolors='black', linewidths=1.0)
       #------------
       networkx.draw_networkx_nodes(graph, pos, nodelist=reactions, node_color='lightcoral', node_size=1000, alpha=1.0, 
                                   node_shape='s',  edgecolors='black', linewidths=1.0)
       #------------
       networkx.draw_networkx_edges(graph, pos, arrows=True, arrowsize=25, width=2)
       #------------
       networkx.draw_networkx_labels(graph, pos)
       #------------
       plt.show(block=False)
    #====================
    return species, reactions, graph
#=================================================================================================================================
#=================================================================================================================================
"""
This Python function derives our new Laplacian model from its stoichiometric model.
"""
#----------------------------------------------------------------------------
def Laplacian_model (stoich, species_sym, rates_sym, sr_graph):
    #-------------------------------
    def incidence_matrix(my_digraph):
        #---
        #The list of nodes and edges
        nodes = list(my_digraph.nodes())
        edges = list(my_digraph.edges())
        #---
        # Create an incidence matrix with shape (number of nodes, number of edges)
        incidence_matrix = numpy.zeros((len(nodes), len(edges)))
        #---
        #Fill the incidence matrix
        for edge_idx, (u, v) in enumerate(edges):
            u_idx = nodes.index(u)
            v_idx = nodes.index(v)
            # Outgoing edge from u: -1, incoming edge to v: +1
            incidence_matrix[u_idx, edge_idx] = -1
            incidence_matrix[v_idx, edge_idx] = 1
        #---
        return incidence_matrix
    #-------------------------------
    Omega = numpy.where(numpy.minimum(numpy.array(stoich), 0) == -0.0, 0.0, - numpy.minimum(numpy.array(stoich), 0))
    #---
    log_species_sym = [sympy.log(item) for item in species_sym]
    expr_sym = numpy.dot(numpy.transpose(Omega), numpy.transpose(log_species_sym))
    exp_expr_sym = [sympy.exp(item) for item in expr_sym]
    #---
    monomial_sym = [sympy.simplify(item) for item in exp_expr_sym]
    #---
    rationals_sym = [sympy.simplify(rates_sym[j]/monomial_sym[j]) for j in range(len(rates_sym))]
    #-------------------------------
    incidence = incidence_matrix(sr_graph)
    #-------------------------------
    weights_1 = [numpy.dot(rationals_sym[j], sympy.diff(monomial_sym[j], species_sym[i]))  for i in range(len(stoich)) for j in range(len(stoich[0])) if stoich[i][j] < 0]
    weights_2 = [stoich[i][j]  for i in range(len(stoich)) for j in range(len(stoich[0])) if stoich[i][j] > 0]
    weights = weights_1 + weights_2
    #-------------------------------
    Gamma = sympy.Matrix(numpy.diag(weights))
    #-------------------------------
    Delta = numpy.minimum(numpy.array(incidence), 0)
    #-------------------------------
    Laplacian = incidence*Gamma*numpy.transpose(Delta)
    #-------------------------------
    d = - numpy.sum(stoich, axis=0)
    #-------------------------------
    r = sympy.Matrix(sympy.symbols(f'r1:{len(rates_sym)+1}'))
    #---                      
    x = sympy.Matrix(species_sym.copy() + list(r))
    #-------------------------------
    Lambda = (-Laplacian * x).tolist()
    #-------------------------------
    return r, d, Laplacian, Lambda
#=================================================================================================================================
#=================================================================================================================================
"""
This Python function derives the stoichiometric model from our new Laplacian model.
"""
#----------------------------------------------------------------------------
def Laplacian_ode (r, d, Laplacian, Lambda, species, constant_species_sym, constant_species_values, parameter_sym, parameter_values, compartment_sym):
    #---
    n_r = len(d)
    #---
    n_s = len(Lambda) - n_r
    #---
    Gamma = sympy.Matrix(numpy.diag(d))
    Lambda_2 = sympy.Matrix([Lambda[j] for j in range(n_s, n_s + n_r)])
    #---
    equations = sympy.Eq(Gamma*r, Lambda_2)
    #---
    Laplacian_rates_ = sympy.solve(equations, r)
    Laplacian_rates = list(Laplacian_rates_.values())
    #---
    L_11 = Laplacian[range(n_s), range(n_s)]
    L_12 = Laplacian[range(n_s), range(n_s,n_s+n_r)]
    L_21 = Laplacian[range(n_s,n_s+n_r), range(n_s)]
    L_22 = Laplacian[range(n_s,n_s+n_r), range(n_s,n_s+n_r)]
    #---
    v = -((L_22+Gamma).inv())*L_21*sympy.Matrix(species)  
    #---
    ds_dt = - L_11*sympy.Matrix(species) - L_12*v
    ds_dt_sym = list(ds_dt)
    ds_dt_sym = [sympy.expand(item) for item in ds_dt_sym]
    #---
    ds_dt = sympy.Matrix([expr.subs({ci: 1 for ci in compartment_sym}) for expr in ds_dt])
    ds_dt =  sympy.Matrix([expr.subs({parameter_sym[i]: parameter_values[i] for i in range(len(parameter_sym))}) for expr in ds_dt])
    ds_dt =  sympy.Matrix([expr.subs({constant_species_sym[i]: constant_species_values[i] for i in range(len(constant_species_sym))}) for expr in ds_dt])
    ds_dt = [sympy.simplify(sympy.expand(item)) for item in ds_dt] 
    #---
    return ds_dt_sym, ds_dt, Laplacian_rates
#=================================================================================================================================
#=================================================================================================================================
"""
This Python function applies Kron reduction to the chemical reaction network described by a species-reaction graph and modeled using our new hybrid system.
"""
#----------------------------------------------------------------------------
def Kron_reduction (species_sym, d, r, Laplacian, Lambda, I_del, J_del):
    #---------------
    I_del.sort() #Species indices for deletion
    I_red = [i for i in range(len(species_sym)) if i not in I_del] #Species indices to keep
    species_sym_del = [species_sym[i] for i in I_del] #Species names for deletion
    species_sym_red = [species_sym[i] for i in I_red] #Species names to keep
    #---------------
    J_del.sort() #Reactions indices for deletion
    J_red = [i for i in range(len(r)) if i not in J_del] #Reactions indices to keep
    reactions_sym_del = [r[i] for i in J_del] #Reactions names for deletion
    reactions_sym_red = [r[i] for i in J_red] #Reactions names to keep
    #---------------
    d_ = [d[i] for i in J_red]
    d_red = numpy.array(d_) 
    x = species_sym_red + reactions_sym_red
    #---------------
    IJ_del = I_del + [j + len(species_sym) for j in J_del] #Node indices to delete
    IJ_red = I_red + [j + len(species_sym) for j in J_red] #Node indices to keep
    #---------------
    s_st_expr = sympy.solve(sympy.Eq(sympy.Matrix([Lambda[i] for i in I_del]), sympy.Matrix([0 for i in I_del])), species_sym_del)
    #---
    if isinstance(s_st_expr, dict):
        s_st = list(s_st_expr.values())[:len(I_del)]
    elif isinstance(s_st_expr, list):
        s_st = s_st_expr
    #---    
    if s_st:
        L = Laplacian.xreplace({species_sym_del[i]: s_st[i] for i in range(len(species_sym_del))})
    else:
        L = Laplacian   
    #---------------
    L_11 = L[IJ_red, IJ_red]
    L_12 = L[IJ_red, IJ_del]
    L_21 = L[IJ_del, IJ_red]
    L_22 = L[IJ_del, IJ_del]
    #---------------
    L_22_inv = L_22.pinv()
    Schur = (L_11 - L_12*L_22_inv*L_21).applyfunc(nsimplify).applyfunc(simplify)   
    #---------------
    Lambda_red = - Schur* sympy.Matrix(x)
    #---------------
    #---------------
    return sympy.Matrix(reactions_sym_red), d_red, Schur, Lambda_red
#=================================================================================================================================
#=================================================================================================================================
"""
This Python function visualizes the comparison of concentration plots between the original and reduced models.
"""
#----------------------------------------------------------------------------
def comparison_plots (time_datapoint, concentrations_original, species_sym, tau, concentrations_reduced, reduced_index, species_compare, file_name, verbose):
    #--------------
    if verbose:
        #------------
        if not species_compare:
               species_compare = [i for i in range(len(concentrations_original))]
        #------------
        for i in range(len(species_compare)):
            #--------------------------------
            fig = plt.figure()
            fig.canvas.manager.set_window_title(f"Concentration plots: {file_name}: {species_sym[species_compare[i]]}")
            #--------------------------------
            tau_i_index = min([j for j in range(len(time_datapoint)-1) if time_datapoint[j] <= tau[i] and time_datapoint[j+1] >= tau[i]])
            #--------------------------------            
            if len(concentrations_original):            
               s_or_i = concentrations_original[i]
               plt.plot(time_datapoint[:tau_i_index], s_or_i[:tau_i_index], label='Original model', color='blue')
            #--------------------------------  
            if len(concentrations_reduced):
                #--------------
                if reduced_index == 0:
                   label_name = f'Reduced model'
                else:
                    label_name = f'Reduced model {reduced_index+1}' 
                #--------------    
                s_red_i = concentrations_reduced[i]                
                plt.plot(time_datapoint[:tau_i_index], s_red_i[:tau_i_index], label=label_name, linestyle='--', color='red')
                #--------------
                #plt.plot([], [], label=f'Error integral = {Error_integral_values[reduced_index]:.2f}', color='white')
            #--------------------------------  
            plt.xlabel('Time', fontsize=16)
            plt.ylabel('Concentration', fontsize=16,)
            plt.xticks(fontsize=12)  # Adjust fontsize for x-axis
            plt.yticks(fontsize=12)
            plt.title(f'{species_sym[species_compare[i]]}', fontsize=20)
            plt.legend(prop={'size': 14})
            plt.show(block=False) 
    #--------------
    return
#=================================================================================================================================
#=================================================================================================================================
"""
This Python function checks whether a given bipartite graph is a valid species-reaction graph. 
"""
#----------------------------------------------------------------------------
def is_species_reaction_graph(species_nodes, reaction_nodes, graph):
        #--------------
        test_1 = graph.number_of_edges() == 0 
        test_2 =  any(graph.has_edge(u, v) for u in species_nodes for v in species_nodes)
        test_3 =  any(graph.has_edge(u, v) for u in reaction_nodes for v in reaction_nodes)
        test_4 = any(graph.has_edge(u, v) and graph.has_edge(v, u) for u in species_nodes for v in reaction_nodes)
        test_5 = any((graph.in_degree(u) + graph.out_degree(u)) == 0 for u in reaction_nodes)
        test_6 = any((graph.in_degree(u) + graph.out_degree(u)) == 0 for u in species_nodes)
        #--------------
        if test_1 or test_2 or test_3 or test_4 or test_5 or test_6:
            return False
        else:
            return True
#=================================================================================================================================
#=================================================================================================================================
"""
This Python function removes a specified species and reaction from the given species-reaction graph.
"""
#----------------------------------------------------------------------------
def remove_species_reaction (species_G, reactions_G, G, species, reaction, verbose ,file_name):
    #--------
    #--------
    G1 = G.copy()
    species_G1 = species_G.copy()
    reactions_G1 = reactions_G.copy()
    #--------
    #--------
    new_edges = []    
    #---
    if species:
        neighborhs = list(G1.successors(species)) + list(G1.predecessors(species))
        #---
        for neigh_1 in neighborhs:
            for neigh_2 in neighborhs:
                if G1.has_edge(neigh_1, species) and G1.has_edge(species, neigh_2) and neigh_1 != neigh_2:
                    new_edges.append((neigh_1, neigh_2)) 
        #---
        G1.remove_node(species)
        #---
        species_G1.remove(species)
    #--------
    if reaction:
        neighborhs = list(G1.successors(reaction)) + list(G1.predecessors(reaction))
        #---
        for neigh_1 in neighborhs:
            for neigh_2 in neighborhs:
                if G1.has_edge(neigh_1, reaction) and G1.has_edge(reaction, neigh_2) and neigh_1 != neigh_2:
                    new_edges.append((neigh_1, neigh_2)) 
        #---
        G1.remove_node(reaction)
        #---
        reactions_G1.remove(reaction) 
    #--------
    #--------
    G1.add_edges_from(new_edges)
    #--------
    if verbose:
        #------------
       fig = plt.figure()
       fig.canvas.manager.set_window_title(f"Species-reaction graph: {file_name}")
       #------------
       pos = networkx.spring_layout(G1)
       #pos = networkx.circular_layout(graph)
       #------------
       networkx.draw_networkx_nodes(G1, pos, nodelist=species_G1, node_color='lightblue', node_size=1000, alpha=1.0, 
                                   node_shape='o',  edgecolors='black', linewidths=1.0)
       #------------
       networkx.draw_networkx_nodes(G1, pos, nodelist=reactions_G1, node_color='lightcoral', node_size=1000, alpha=1.0, 
                                   node_shape='s',  edgecolors='black', linewidths=1.0)
       #------------
       networkx.draw_networkx_edges(G1, pos, arrows=True, arrowsize=25, width=2)
       #------------
       networkx.draw_networkx_labels(G1, pos)
       #------------
       plt.show(block=False)
    #====================

    #--------
    return species_G1, reactions_G1, G1   
#=======================================================================================
#=======================================================================================
#=================MODEL REDUCTION=======================================================
#=======================================================================================
if __name__ == '__main__':
    #--------------------------------------------
    file_path, file_name = get_file()
    sbml_model, sbml_document = get_model(file_path)
    #----------
    species_sym, initial_conc, const_species_sym, const_species_values = model_species (sbml_model)
    #----------
    stoich_matrix = stoichiometric_matrix (sbml_model)   
    #----------
    parameter_sym, parameter_values = model_parameters(sbml_model)
    #----------
    compartment_sym, compartment_size = model_compartments(sbml_model)
    #----------
    rates_sym = model_rates (sbml_model)
    #-----------------------------------------------------------
    verbose_1 = 0
    sr_species, sr_reactions, sr_graph = species_reaction_graph(stoich_matrix, species_sym, file_name, verbose_1)
    #-----------------------------------------------------------
    ds_dt_sym, ds_dt = stoichiometric_ode_model (stoich_matrix, const_species_sym, const_species_values, rates_sym, parameter_sym, parameter_values, compartment_sym)
    #-----------------------------------------------------------
    T = 180
    n = 1000
    delta = 0.1
    method_int = 'LSODA'
    #---
    steady_state, tau = equilibrium (ds_dt, T, n, species_sym, initial_conc, delta, method_int)
    #-----------------------------------------------------------
    verbose_2 = 0
    stoich_time, stoich_concentrations = ode_model_integrations (ds_dt, max(tau), n, species_sym, initial_conc, method_int, file_name, verbose_2)
    #-----------------------------------------------------------
    start_time = time.time()
    #-----------------------------------------------------------
    r, d, Laplacian, Lambda = Laplacian_model (stoich_matrix, species_sym, rates_sym, sr_graph)
    #-----------------------------------------------------------
    ds_dt_lp_sym, ds_dt_lp, Laplacian_rates = Laplacian_ode (r, d, Laplacian, Lambda, species_sym, const_species_sym, const_species_values, parameter_sym, parameter_values, compartment_sym)
    #-----------------------------------------------------------
    #I_del = [1, 2, 5, 7, 10]
    #J_del = [0, 1, 4, 5, 7, 8, 12]
    I_del = [8, 10, 11, 12, 13]
    J_del = [4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 16, 17]
    #---
    r_red, d_red, Schur, Lambda_red = Kron_reduction (species_sym, d, r, Laplacian, Lambda, I_del, J_del)
    species_red_sym = [species_sym[i] for i in range(len(species_sym)) if i not in I_del]
    #-----------------------------------------------------------
    #-----------------------------------------------------------
    ds_dt_sym_red, ds_dt_red, Laplacian_rates_red = Laplacian_ode (r_red, d_red, Schur, Lambda_red, species_red_sym, const_species_sym, const_species_values, parameter_sym, parameter_values, compartment_sym)
    #-----------------------------------------------------------
    verbose_3 = 0
    initial_conc_red = [initial_conc[i] for i in range(len(initial_conc)) if i not in I_del]
    red_time, red_concentrations = ode_model_integrations (ds_dt_red, max(tau), n, species_red_sym, initial_conc_red, method_int, file_name, verbose_3)
    #-----------------------------------------------------------
    I_red = [i for i in range(len(initial_conc)) if i not in I_del]
    key_species = [0,3,4]
    key_species = I_red
    key_species_red = [item for item in I_red if item in key_species]
    key_species_red_index = [I_red.index(item) for item in  key_species_red]
    #---
    x = stoich_concentrations[key_species,:]
    #---
    y = red_concentrations[key_species_red_index,:]
    #--------------- 
    error_value = error_integral_value = simmetrized_error_integral (x, y, tau)
    #-----------------------------------------------------------
    tau_ = [200 for _ in tau]
    stoich_time_2, stoich_concentrations_2 = ode_model_integrations (ds_dt, max(tau_), n, species_sym, initial_conc, method_int, file_name, verbose_2)
    red_time_2, red_concentrations_2 = ode_model_integrations (ds_dt_red, max(tau_), n, species_red_sym, initial_conc_red, method_int, file_name, verbose_3)
    #---
    x_2 = stoich_concentrations[key_species,:]
    #---
    y_2 = red_concentrations[key_species_red_index,:]
    #--------------- 
    verbose_4 = 1
    comparison_plots (stoich_time_2, x_2, species_red_sym, tau_, y_2, 0, [], file_name, verbose_4)
    #---------------
    species_G = sr_species
    reactions_G = sr_reactions 
    G = sr_graph
    #---
    for i in I_del:
        species_ = sr_species[i]
        species_G, reactions_G, G = remove_species_reaction (species_G, reactions_G, G, species_, [], 0 ,file_name)
    #---    
    for i in J_del:
        reaction__ = sr_reactions[i]
        species_G, reactions_G, G = remove_species_reaction (species_G, reactions_G, G, [], reaction__, 0 ,file_name)
    #-----------------------------------------------------------
    #-----------------------------------------------------------
    end_time = time.time()
    execution_time = end_time - start_time
    #-----------------------------------------------------------
    #-----------------------------------------------------------
    print()
    print()   
    print(f"{BOLD}{RED}==============================================={RESET}{RESET}")
    print(f"{BOLD}{RED}==============================================={RESET}{RESET}")
    print(f"{BOLD}{RED}=== Initiating Data Output Handling Process ==={RESET}{RESET}")
    print()
    print() 
    print(f"{BOLD} BioModel ID: {RESET}{BLUE}{file_name}{RESET}")
    print()
    print(f"{BOLD}\t- List of species:{RESET} {BLUE}{species_sym}{RESET}")
    print(f"{BOLD}\t- List of species initial values:{RESET} {BLUE}{initial_conc}{RESET}")
    print()
    if const_species_sym:
       print(f"{BOLD}\t- List of constant species:{RESET} {BLUE}{const_species_sym}{RESET}")
       print(f"{BOLD}\t- List of constant species values:{RESET} {BLUE}{const_species_values}{RESET}")
       print()
    print(f"{BOLD}\tStoichiometric matrix: {RESET}{BLUE}{stoich_matrix}{RESET}")
    print()
    print(f"{BOLD}\t- List of parameters:{RESET} {BLUE}{parameter_sym}{RESET}")
    print(f"{BOLD}\t- List of parameter valuess:{RESET} {BLUE}{parameter_values}{RESET}")
    print()
    print(f"{BOLD}\t- List of compartments:{RESET} {BLUE}{compartment_sym}{RESET}")
    print(f"{BOLD}\t- List of compartment sizes:{RESET} {BLUE}{compartment_size}{RESET}")
    print()
    print(f"{BOLD}\t- List of reaction rates:{RESET} {BLUE}{rates_sym}{RESET}")
    print()
    print(f"{BOLD}\t- List of edges in the SR graph:{RESET} {BLUE}{sr_graph.edges}{RESET}")
    print()
    print(f"{BOLD}\t- List of stoichiometric ODEs:{RESET} {BLUE}{ds_dt_sym}{RESET}")
    print()
    print(f"{BOLD}\t- List of equilibrium points: {RESET} {BLUE}{steady_state}{RESET}")
    print(f"{BOLD}\t- List of species equilibrium times: {RESET} {BLUE}{tau}{RESET}")
    print(f"{BOLD}\t- Equilibrium time: {RESET} {BLUE}{max(tau)}{RESET}")
    print()
    print(f"{BOLD}\t- List of symbolic rates: {RESET} {BLUE}{r.tolist()}{RESET}")
    print(f"{BOLD}\t- List of net stoichiometric imbalancees: {RESET} {BLUE}{list(d)}{RESET}")
    print(f"{BOLD}\t- Laplacian matrix: {RESET} {BLUE}{Laplacian.tolist()}{RESET}")
    print(f"{BOLD}\t- Laplacian model: {RESET} {BLUE}{Lambda}{RESET}")
    print()
    print(f"{BOLD}\t- List of Laplacian ODEs: {RESET} {BLUE}{ds_dt_lp_sym}{RESET}")
    print(f"{BOLD}\t- List of Laplacian rates: {RESET} {BLUE}{Laplacian_rates}{RESET}")
    print()
    print(f"{BOLD}\t- List of reduced Laplacian ODEs: {RESET} {BLUE}{ds_dt_sym_red}{RESET}")
    print(f"{BOLD}\t- List of reduced Laplacian rates: {RESET} {BLUE}{Laplacian_rates_red}{RESET}")
    print()
    print(f"{BOLD}\t- The value of the symmetrized error integral: {RESET} {BLUE}{error_value}{RESET}")
    print()
    if is_species_reaction_graph(species_G, reactions_G, G):
       print(f"{BOLD}\t- The resulting graph is a vaid species-reaction graph:{RESET}")
       print(f"{BOLD}\t- List of edges in the SR graph:{RESET} {BLUE}{G.edges}{RESET}")
       print() 
    else:
       print(f"{BOLD}\t- The resulting graph is a NOT vaid species-reaction graph:{RESET}")
       print(f"{BOLD}\t- List of edges in the SR graph:{RESET} {BLUE}{G.edges}{RESET}")
       print() 
    print(f"{BOLD}\t- Execution time:{RESET} {BLUE}{execution_time} seconds{RESET}")
    #======================================================================================= 
    #==============SAVING THE OUTPUT======================================================== 
    selected_outputs = {
        "Species ID": str(species_sym),
        "Buffered Species ID": const_species_sym,
        "Initial Concentrations": initial_conc,
        "Constant Concentrations": const_species_sym,
        "Parameter IDs": str(parameter_sym),
        "Parameter Values": parameter_values,
        "Compartment ID": str(compartment_sym),
        "Compartment Size": compartment_size,
        "Stoichiometric matrix": str(stoich_matrix),
        "Original model": str(ds_dt_sym),
        "Removed species": str([item for item in species_sym if item not in species_red_sym]),
        "Remaining species": str(species_red_sym),
        "Reduced rates": str(Laplacian_rates_red),
        "Reduced model": str(ds_dt_sym_red),
        "Value of the symmertic error integral": error_integral_value,
    }
    #------------------------------------------------------------------------------
    current_dir = os.path.dirname(os.path.abspath(__file__))
    json_file_path = os.path.join(current_dir, str(file_name) + ".json")
    #Save the outputs to the json file
    with open(json_file_path, "w") as json_file:
        json.dump(selected_outputs, json_file, indent=4)
        json_file.write('\n')
    print(f"{BOLD}\t- The outputs successfully saved to::{RESET} {BLUE}{json_file_path} seconds{RESET}")
    #=================================================================================================================================
    #=================================================================================================================================
    if verbose_1 or verbose_2 or verbose_3 or verbose_4:
       input()
