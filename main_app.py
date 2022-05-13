from src.reduction import ReductionMethod
import cantera as ct


def main() -> None:
    """
    Edit all the variables in this function to perform the reduction.

    Right now, it has DRG and DRGEP. Put all the state files in a folder and pass
    the folder path to load the condition.

    The automation of the result mechansim and testing will be implemented later.

    Returns a set containing all the species deemed important.
    """

    # detailed mechanism
    detailed_gas = ct.Solution("gri30.xml")

    # folder path (use \\ to separate the folders)
    file_dir = "folder_path\\"

    # Reduction parameters
    threshold = 0.05
    important_species = [ "OH", "CH4", "O2", "N2", "CO2", "H2O", "CO"]
    reduction_method = "DRGEP"            # accepts "DRG" or "DRGEP"

    # call the reduction
    red = ReductionMethod(detailed_gas, det_dir=file_dir)
    final_spc_list, rii_list = red.run_reduction(reduction_method, threshold, important_species)
    
    print(f" Final red_mech contains {len(final_spc_list)} species.\n (red/det) = ({len(final_spc_list)}/{detailed_gas.n_species})\n")
    for spc in final_spc_list:
        print(spc)





if __name__ == "__main__":
    main()

