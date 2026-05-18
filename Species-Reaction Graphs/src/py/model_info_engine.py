import libsbml, sys, numpy, networkx, re, sympy
from sympy import sympify
from libsbml import *
from scipy.integrate import solve_ivp
from tkinter import filedialog, Tk
import matplotlib.pyplot as plt
#=======================================================================================================================
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"           # Reaction ID
#-----------------------------------------------------------------
numpy.set_printoptions(threshold=numpy.inf) # Display the full array
#=======================================================================================================================
def get_file():
    root = Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        parent=root,
        title="Select SBML file",
        filetypes=[("All files", "*.*")]
    )

    if not file_path:
        print(f"\n{BOLD}No file selected.\n")
        sys.exit(1)

    file_name = os.path.basename(file_path).split('.')[0]
    return file_path, file_name
#=======================================================================================================================
def get_model(file_path: str) -> tuple[libsbml.Model, libsbml.SBMLDocument]:
    #---
    reader = libsbml.SBMLReader()
    #---
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
#=======================================================================================================================
def get_model_species (model):
    #---
    buffered_species_id = [species.getId() for species in model.getListOfSpecies() if species.getConstant()]
    bufferred_values = [species.getInitialAmount() or species.getInitialConcentration() for species in model.getListOfSpecies() if species.getId() in buffered_species_id]
    #---
    species_id = [species.getId() for species in model.getListOfSpecies() if not species.getConstant()]
    init_conc = [species.getInitialAmount() or species.getInitialConcentration() for species in model.getListOfSpecies() if species.getId() in species_id]
    #---
    return buffered_species_id, bufferred_values, species_id, init_conc
#=======================================================================================================================
def get_model_reactions (model):
    #---
    #The reaction IDs
    reactions_id = [reaction.getId() for reaction in model.getListOfReactions()]
    #---
    reactions_info_species_id = [
            f"{' + '.join([f'{s.getStoichiometry()}*{model.getSpecies(s.getSpecies()).getId()}' for s in reaction.getListOfReactants()])} "
            f"{'<--->' if reaction.getReversible() else '--->'} "
            f"{' + '.join([f'{s.getStoichiometry()}*{model.getSpecies(s.getSpecies()).getId()}' for s in reaction.getListOfProducts()])}"
            for reaction in model.getListOfReactions()
            ]
    #---
    clean_reactions_info_species_id = [re.sub(r'(\d+)\.0\*([A-Za-z]+)', lambda m: f'{m.group(2)}' if m.group(1) == '1' else f'{m.group(1)}*{m.group(2)}', reaction) for reaction in reactions_info_species_id]
    #---
    reactions_info_irrev = []    
    reactions_id_irrev = []
    #--------------------------------------------             
    for reaction, reaction_id in zip(clean_reactions_info_species_id, reactions_id):    
        #-----------------------    
        if '<--->' in reaction:
                #---
                substrate = reaction.split('<--->')[0].strip()
                product = reaction.split('<--->')[1].strip()
                #---
                reaction_if = f"{substrate} ---> {product}"
                reaction_if_id = reaction_id + "_f"
                #---
                reaction_ir = f"{product} ---> {substrate}"
                reaction_ir_id = reaction_id + "_r"
                #---
                reactions_info_irrev.append(reaction_if)
                reactions_info_irrev.append(reaction_ir)
                #---       
                reactions_id_irrev.append(reaction_if_id)
                reactions_id_irrev.append(reaction_ir_id)
        else:
                reactions_info_irrev.append(reaction)  
                reactions_id_irrev.append(reaction_id)
        #-----------------------
    return reactions_id_irrev, reactions_info_irrev
#=======================================================================================================================
def get_model_rates (model):

    #---
    rates = [reaction.getKineticLaw().getFormula() if reaction.getKineticLaw() else "No rate law" 
             for reaction in model.getListOfReactions()]
    #---
    rates_irrev = []
    #------------------------
    for rate in rates:
        if '-' in rate:    
           rate_f = rate.split('-')[0].strip()
           rate_r = rate.split('-')[1].strip()
           rates_irrev.append(rate_f)
           rates_irrev.append(rate_r)
        else:
           rates_irrev.append(rate) 
    #------------------------
    return rates_irrev
#=======================================================================================================================
def get_model_parameters (model):
    #---
    parameter_id = [parameter.getId() for parameter in model.getListOfParameters()]
    parameter_values = [parameter.getValue() for parameter in model.getListOfParameters()]
    #---
    return parameter_id, parameter_values
#=======================================================================================================================
def get_model_compartments (model):
    #---
    compartment_id = [compartment.getId() for compartment in model.getListOfCompartments()]
    compartment_size = [compartment.getSize() for compartment in model.getListOfCompartments()]
    #---
    return compartment_id, compartment_size
#=======================================================================================================================
def get_stoichiometric_matrix (species_id, reactions):
    #----------------------------------------------------------------------------
    reactions_sym = []
    for reaction in reactions:
        #---------------------------------------------------
        if '--->' in reaction:
            #---
            substrate = reaction.split('--->')[0].strip()
            product = reaction.split('--->')[1].strip()
        #---------------------------------------------------
        if substrate and product:
           reaction_expression = sympify(product) - sympify(substrate)
           reactions_sym.append(reaction_expression)
        #---------------------------------------------------
        if substrate and not product:
           reaction_expression = - sympify(substrate)
           reactions_sym.append(reaction_expression) 
        #---------------------------------------------------
        if not substrate and product:
           reaction_expression = sympify(product) 
           reactions_sym.append(reaction_expression)
    #----------------------------------------------------------------------------
    species_sym = sympify(species_id)
    stoich = [[reaction_sym.coeff(species_sym_i) for reaction_sym in reactions_sym] for species_sym_i in species_sym]
    #----------------------------------------------------------------------------
    return stoich
#=======================================================================================================================
def sr_graph(stoichiometric, species, file_name, verbose=0):
    #------------
    graph = networkx.DiGraph()
    #------------
    #Label the species and the reaction nodes
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
    if verbose:
       #------------
       fig = plt.figure()
       fig.canvas.manager.set_window_title(f"Species-reaction graph: {file_name}")
       #------------
       #pos = networkx.spring_layout(graph)
       pos = networkx.circular_layout(graph)
       #pos = graphviz_layout(graph, prog="dot")
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
#=======================================================================================================================
def ode_model (stoich, bufferred_sym, bufferred_values, rates_sym, parameter_sym, parameter_values, compartment_sym):
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    S = sympy.Matrix(stoich)
    #---
    v1 = sympy.Matrix([expr.subs({ci: 1 for ci in compartment_sym}) for expr in rates_sym.copy()])
    #---
    ds_dt_sym = list(S*v1)
    #---
    v2 =  sympy.Matrix([expr.subs({parameter_sym[i]: parameter_values[i] for i in range(len(parameter_sym))}) for expr in v1.copy()])
    #---
    v =  sympy.Matrix([expr.subs({bufferred_sym[i]: bufferred_values[i] for i in range(len(bufferred_sym))}) for expr in v2.copy()])
    #---
    ds_dt = S*v
    #------------------------------------------------------------------
    #------------------------------------------------------------------  
    return ds_dt_sym, ds_dt
#=======================================================================================================================
def ode_integrations (ds_dt, T, n, species_sym, init_conc, method_int, file_name, verbose):
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
    solutions_final = solve_ivp(fun_int,  (0, T), init_conc, t_eval=numpy.linspace(0, T, n), method=method_int)
    solution_concentrations = solutions_final.y
    solution_time = solutions_final.t
    #------------------------------------------------------------------
    if verbose:
        for i in range(len(species_sym)):
            # Create a new figure for each plot
            fig = plt.figure()
            fig.canvas.manager.set_window_title(f"Concentration plot: {file_name}")  # Set window title to species name
            
            # Plot the solution for the current species
            plt.plot(solution_time, solution_concentrations[i], label=f'{species_sym[i]}')
            plt.xlabel('Time')
            plt.ylabel('Concentration')
            plt.title(f'{species_sym[i]}')
            plt.legend()
            plt.show(block=False) 
    #------------------------------------------------------------------         
    return solution_time, solution_concentrations 
#=======================================================================================================================
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
#=======================================================================================================================
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
    v = -((L_22+Gamma).pinv())*L_21*sympy.Matrix(species)  
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
#=======================================================================================================================
