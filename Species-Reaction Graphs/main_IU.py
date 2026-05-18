import numpy, sympy
from sympy import sympify
#---
from model_info_engine import  (
                    get_model, 
                    get_file,
                    get_model_species,
                    get_model_reactions,
                    get_model_rates,
                    get_model_parameters,
                    get_model_compartments,
                    get_stoichiometric_matrix,
                    sr_graph,
                    ode_model,
                    ode_integrations,
                    Laplacian_model,
                    Laplacian_ode
                    )    
#------------------------------------------------------------------------------------------------------------------------
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"           # Reaction ID
#------------------------------------------------------------------------------------------------------------------------
numpy.set_printoptions(threshold=numpy.inf) # Display the full array
#=======================================================================================================================
if __name__ == '__main__':
    ################################################################
    ################################################################
    #Get files paths and files names
    file_path, file_name = get_file()
    print(f"{BOLD}\n\n\n\tSelected BioModel: {RESET}{BLUE}{file_name}{RESET}\n")
    #===========
    #Model and document
    sbml_model, sbml_document = get_model(file_path)   
    #===========
    #Constant species and Species
    buffered_species_id, bufferred_values, species_id, init_conc = get_model_species (sbml_model)
    print(f"{BOLD}\tBuffered species: {RESET}{BLUE}{buffered_species_id}{RESET}\n")    
    print(f"{BOLD}\tBuffered species values: {RESET}{BLUE}{bufferred_values}{RESET}\n")   
    print(f"{BOLD}\tSpecies: {RESET}{BLUE}{species_id}{RESET}\n") 
    print(f"{BOLD}\tInitial concentrations: {RESET}{BLUE}{init_conc}{RESET}\n")  
    #=========== 
    #Reactions in terms of species
    reactions_id_irrev, reactions_info_irrev = get_model_reactions (sbml_model)
    print(f"{BOLD}\tReaction IDs: {RESET}{BLUE}{reactions_id_irrev}{RESET}\n")  
    print(f"{BOLD}\tReactions: {RESET}{BLUE}{reactions_info_irrev}{RESET}\n")  
    #=========== 
    #Reaction rates
    rates_irrev = get_model_rates (sbml_model)
    print(f"{BOLD}\tReaction rates: {RESET}{BLUE}{rates_irrev}{RESET}\n")  
    #=========== 
    #Parameters and its values
    parameter_id, parameter_values = get_model_parameters (sbml_model)
    print(f"{BOLD}\tParameters: {RESET}{BLUE}{parameter_id}{RESET}\n")  
    print(f"{BOLD}\tParameter values: {RESET}{BLUE}{parameter_values}{RESET}\n")  
    #=========== 
    #Compartments
    compartment_id, compartment_size = get_model_compartments (sbml_model)
    print(f"{BOLD}\tCompartments: {RESET}{BLUE}{compartment_id}{RESET}\n")  
    print(f"{BOLD}\tCompartment sizes: {RESET}{BLUE}{compartment_size}{RESET}\n")  
    #=========== 
    #Stoichiometric matrix
    stoich = get_stoichiometric_matrix (species_id, reactions_info_irrev)
    print(f"{BOLD}\tStoichiometric matrix: {RESET}{BLUE}{numpy.array(stoich)}{RESET}\n")  
    #=========== 
    #Species reaction graph
    species, reactions, graph = sr_graph(stoich, species_id, file_name, verbose=0)
    #=========== 
    species_sym = sympify(species_id)
    bufferred_sym = sympify(buffered_species_id)
    rates_sym = sympify(rates_irrev)
    parameter_sym = sympify(parameter_id)
    compartment_sym =sympify(compartment_id)
    #=========== 
    #ODE model
    ds_dt_sym, ds_dt = ode_model (stoich, bufferred_sym, bufferred_values, rates_sym, parameter_sym, parameter_values, compartment_sym)
    print(f"{BOLD}\tODE model: {RESET}{BLUE}{numpy.array(ds_dt)}{RESET}\n")  
    #===========
    #Model simulation
    T = 100000
    n = 100000
    method_int = 'LSODA'
    solution_time, solution_concentrations = ode_integrations (ds_dt, T, n, species_sym, init_conc, method_int, file_name+"_original", verbose=0)
    #===========
    #Equilibrium
    tau = [10^5, 10^4, 2*10^3, 3*10^3, 4*10^3, 2*10^3, 2*10^3, 500, 500, 6*10^4, 500, 500, 500, 500]
    #===========
    r, d, Laplacian, Lambda = Laplacian_model (stoich, species_sym, rates_sym, graph)
    #===========
    ds_dt_L_sym, ds_dt_L, Laplacian_rates = Laplacian_ode (r, d, Laplacian, Lambda, species, bufferred_sym, bufferred_values, parameter_sym, parameter_values, compartment_sym)
    print(f"{BOLD}\tLaplacian ODE model: {RESET}{BLUE}{numpy.array(ds_dt_L_sym)}{RESET}\n")  
    print(f"{BOLD}\tLaplacian rates: {RESET}{BLUE}{numpy.array(Laplacian_rates)}{RESET}\n") 
    ################################################################
    ################################################################




    ################################################################
    ################################################################
    


            



    # input()








