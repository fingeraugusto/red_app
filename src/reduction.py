from typing import List, Tuple, Set
from .mech_graph import MechGraph
from .utils import get_state
import cantera as ct


class ReductionMethod:
    def __init__(self, det_gas: ct.Solution,  det_dir: str = "", target_conditions_file: str = "",) -> None:
        self.gas = det_gas
        self.threshold = 0.1
        self.all_files_dir = det_dir
        self.target_cond = target_conditions_file
    
    def run_reduction(self, red_method: str, threshold: float, imp_spc_list: List[str]) -> Tuple[Set[str], List[Set[str]]]:
        graph = MechGraph()
        graph.set_cantera_model(self.gas)
        if red_method == "DRG":
            red_function = graph.drg
        else:
            red_function = graph.drgep

        ind_red_mechs = []
        final_red = set()

        i = 0

        # call the reduction in every condition
        for state in get_state(self.all_files_dir):
            print(i)
            i+= 1
            # the 4 element in state is if the file is mass fraction basd
            if state[3]:
                self.gas.TPY = state[0], state[1], state[2]
            else:
                self.gas.TPX = state[0], state[1], state[2]
            
            ind_red_mechs.append(red_function(self.gas, imp_spc_list, threshold))
            final_red |= ind_red_mechs[-1]

            if len(final_red) >= self.gas.n_species:
                print(f"WARNING: no reduction. Aborting the process. Consider chaning the paramaters.")
                return 
        
        return final_red, ind_red_mechs



if __name__ == "__main__":
    gas = ct.Solution("POLIMI_TOT_NOX_1412.CKI.yaml")
    print(gas.n_species)
    # a = input()
    dir = "C:\\Users\\1511 IRON\\Desktop\\Pós\\Doutorado\\PapersToBe\\Virtual_N_HEP\\red_teste2\\"
    red = ReductionMethod(gas, det_dir=dir)
    f, r = red.run_reduction("DRGEP", 0.05, ["CO2", "OH"])
    print(f)
    print(len(f))
