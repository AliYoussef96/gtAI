from math import sin, cos, pi
from gaft import GAEngine
from gaft.components import BinaryIndividual
from gaft.components import Population
from gaft.operators import RouletteWheelSelection
from gaft.operators import TournamentSelection
from gaft.operators import UniformCrossover
from gaft.operators import FlipBitBigMutation

# Built-in best fitness analysis.
from gaft.analysis.fitness_store import FitnessStore
from gaft.analysis.console_output import ConsoleOutput

def gene_algo_corr(dict_tGCN,genetic_code_number,RSCU_df,bacteria, size_pop,generation_number):
  
    # Define population.
    indv_template = BinaryIndividual(ranges=[(0, 1),(0,1),(0,1),(0,1),(0,1)], eps=0.001)
    population = Population(indv_template=indv_template, size=size_pop).init()

    # Create genetic operators.
    selection = RouletteWheelSelection()
    crossover = UniformCrossover(pc=0.5, pe=0.25)
    mutation = FlipBitBigMutation(pm=0.1, pbm=0.5, alpha=0.51)

    # Create genetic algorithm engine.
    # Here we pass all built-in analysis to engine constructor.
    engine = GAEngine(population=population, selection=selection,
                    crossover=crossover, mutation=mutation,
                    analysis=[ConsoleOutput, FitnessStore])


    #########################################
    # Define fitness function.
    @engine.fitness_register
    #@engine.minimize
    def fitness(indv):
        from gtAI import new_flow 

        xug, xci, xai, xgu, xal = indv.solution

        Sug = xug
        Sci = xci
        Sai = xai
        Sgu = xgu
        Sal = xal
        
        Wi_df = new_flow.wi_tai_calc(dict_tGCN ,Sug=xug, Sci=xci, Sai=xai, Sgu=xgu, Sal=xal,genetic_code_number=genetic_code_number ,bacteria=bacteria)
        
        RSCU_wai_corr = new_flow.corr_result(Wi_df,RSCU_df)
        
        corr = float(RSCU_wai_corr["Wi"][0]) 
        if corr:
            return corr
        elif not corr:
            return 0
    
       

    engine.run(ng=generation_number)




