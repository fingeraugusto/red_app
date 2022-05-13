from typing import Dict, List, Set
import cantera as ct
import collections

      
class MechGraph:
    def __init__(self) -> None:
        self.spc_rii = {}
        self.spc_total = {}
        self.spc_to_spc = {}

    def set_cantera_model(self, gas: ct.Solution) -> None:
        species_list = gas.species()
        self.spc_rii = {s.name:0.0 for s in species_list}
        self.spc_total = {s.name:0.0 for s in species_list}
        for i in range(gas.n_reactions):
            self.process_rxn(gas.reaction(i))
            
    def process_rxn(self, rxn: ct.Reaction, rate_value: float = 0.0) -> None:
        for s1 in rxn.reactants:
            self.spc_total[s1] += abs(rxn.reactants[s1] * rate_value)
            if s1 not in self.spc_to_spc:
                self.spc_to_spc[s1] = {}
            for s2 in rxn.products:
                self.spc_total[s2] += abs(rxn.products[s2] * rate_value)
                
                if s2 not in self.spc_to_spc:
                    self.spc_to_spc[s2] = {}
                if s2 not in self.spc_to_spc[s1]:
                    self.spc_to_spc[s1][s2] = 0.0
                if s1 not in self.spc_to_spc[s2]:
                    self.spc_to_spc[s2][s1] = 0.0
                
                self.spc_to_spc[s1][s2] += abs(rxn.reactants[s1] * rate_value)
                self.spc_to_spc[s2][s1] += abs(rxn.products[s2] * rate_value)     
    
    def normalize_rates(self) -> None:
        for s1 in self.spc_to_spc:
            for s2 in self.spc_to_spc[s1]:
                self.spc_to_spc[s1][s2] = self.spc_to_spc[s1][s2] / self.spc_total[s1]

    def reset_spc_values(self) -> None:
        for s in self.spc_rii:
            self.spc_rii[s] = 0.0
            self.spc_total[s] = 0.0
        
        for s1 in self.spc_to_spc:
            for s2 in self.spc_to_spc[s1]:
                self.spc_to_spc[s1][s2] = 0.0
   
    def update_state(self, gas: ct.Solution) -> None:
        self.reset_spc_values()
        rate_of_progress = gas.forward_rate_constants
        for i in range(gas.n_reactions):
            self.process_rxn(gas.reaction(i), rate_value=rate_of_progress[i])
        self.normalize_rates()
    
    def _init_search_queue(self, imp_list: List[str]) -> collections.deque:
        return collections.deque([spc for spc in imp_list if spc in self.spc_to_spc])

    def drg_search(self, imp_list: List[str], threshold: float) -> Set[str]:
        queue = self._init_search_queue(imp_list)
        out_list = set(imp_list)
        for s in imp_list:
            self.spc_rii[s] = 1.0
        
        while queue:
            node_to_proccess = queue.popleft()

            for spc, ri in self.spc_to_spc[node_to_proccess].items():
                self.spc_rii[spc] = max(self.spc_rii[spc], ri)
                if (self.spc_rii[spc] >= threshold) and (spc not in out_list):
                    queue.append(spc)
                    out_list.add(spc)
        
        return out_list  

    def drg(self, gas: ct.Solution, imp_list: List[str], threshold: float) -> Set[str]:
        self.update_state(gas)
        return self.drg_search(imp_list, threshold)

    def drgep_search(self, imp_list: List[str], threshold: float) -> Set[str]:
        queue = self._init_search_queue(imp_list)
        out_list = set(imp_list)
        for s in imp_list:
            self.spc_rii[s] = 1.0
        
        while queue:
            node_to_proccess = queue.popleft()
            out_list.add(node_to_proccess)

            for spc, ri in self.spc_to_spc[node_to_proccess].items():
                _aux = ri * self.spc_rii[node_to_proccess]
                if (spc in out_list) and (_aux > self.spc_rii[spc]):
                    self.spc_rii[spc] = _aux
                    queue.append(spc)
                elif (spc not in out_list):
                    self.spc_rii[spc] = max(self.spc_rii[spc], _aux)
                    if self.spc_rii[spc] >= threshold:
                        queue.append(spc)
                        out_list.add(spc)
        
        return out_list
    
    def drgep(self, gas: ct.Solution, imp_list: List[str], threshold: float) -> Set[str]:
        self.update_state(gas)
        return self.drgep_search(imp_list, threshold)






if __name__ == "__main__":
    gas = ct.Solution("gri30.xml")
    a = MechGraph()
    a.set_cantera_model(gas)
    a.update_state(gas)
    out = {}
    r = a.drg_search(["H2"], 0.1)
    print(r)
    # for s in a.spc_rii:
    #     print(f"{s} - {a.spc_rii[s]}")

    
    