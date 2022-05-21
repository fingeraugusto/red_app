import cantera as ct
from src.reactors import IDTReactor, LFSReactor


def main() -> None:
    """
    Here you will see an implementation of IDT and LFS for a list of conditions.
    The current script wll run and save each condition to the folder, using a easy name convention:
    <model>_<Temperature>_<pressure>_<phi>.csv

    All the species values are molar fractions.
    """
    # set the basic cantera solution
    gas = ct.Solution("gri30.xml")
    fuel = "CH4"
    oxdyzer = "O2:1.0,N2:3.76"
    dir_to_save = "directory_to_save\\"

    # set all T, P, PHI conditions for each type of model. If you do not want to run one type of model, just
    # leave the list empty. It will start with the IDT whch is faster and will move to LFS.
    idt_conditions = [  
        [1000.0, ct.one_atm, 1.0],
        [1100.0, ct.one_atm, 1.0],
        [1200.0, ct.one_atm, 1.0],
        [1300.0, ct.one_atm, 1.0],
        [1400.0, ct.one_atm, 1.0],
        [1500.0, ct.one_atm, 1.0]
                    ]
    
    lfs_conditions = [  
        [350.0, ct.one_atm, 0.8],
        [350.0, ct.one_atm, 0.9],
        [350.0, ct.one_atm, 1.0],
        [350.0, ct.one_atm, 1.1],
        [350.0, ct.one_atm, 1.2],
                    ]
    
    print("\nINFO: Some models could take a long time to run.\n")

    # for every conditon in the IDT list and run the model:
    for condition in idt_conditions:
        T, P, phi = condition
        print(f"INFO: Solving IDT reactor for condition: T = {T:.1f}K, P = {P/ct.one_atm:.1f} atm, phi = {phi:.1f}")

        gas.set_equivalence_ratio(phi, fuel, oxdyzer)
        gas.TP = T, P

        idt_reactor = IDTReactor()
        idt, states, header = idt_reactor.solve(gas, path_to_save=dir_to_save, phi_for_file_name=phi)
    
    # for every conditon in the LFS list and run the model:
    gas.transport_model = "Mix"
    for condition in lfs_conditions:
        T, P, phi = condition
        print(f"INFO: Solving LFS reactor for condition: T = {T:.1f}K, P = {P/ct.one_atm:.1f} atm, phi = {phi:.1f}")

        gas.set_equivalence_ratio(phi, fuel, oxdyzer)
        gas.TP = T, P

        lfs_reactor = LFSReactor()
        lfs, states, header = lfs_reactor.solve(gas, path_to_save=dir_to_save, phi_for_file_name=phi)
    
    print("INFO: All cases solved!!")


#---------------------------------------
if __name__ == "__main__":
    main()