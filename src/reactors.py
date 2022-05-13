import cantera as ct
import numpy as np
from typing import List, Tuple
from scipy import integrate
from copy import copy

"""
Present a simple implementation of IDT reactors and the
cantera implementation of a LFS reactor.


Each model can be called as:

IDT, all_conditions = idt_reactor.solve(gas, flag='T', temp_rise=400)   
IDT, all_conditions = idt_reactor.solve(gas, path_to_save=dir_to_save, phi_for_file_name=phi_value)
WARNINGS: idt_reactor uses only the Temperarure rise. The species peak is still under developement.

LFS, all_conditions = lfs_reactor.solve(gas)
LFS, all_conditions = lfs_reactor.solve(gas, path_to_save=dir_to_save, phi_for_file_name=phi_value)

The saved file will be named as:
<model_type>_<TEMP in K>_<PRESSURE in atm>_<PHI, default is 1.0>.csv

The first line will contain the propertie value (IDT or LFS).
The second line will be a header containing specifics and the conditions.
The rest of the file will present the values in a csv format.
"""

# --------------------------------------------------------------------------------------------------
# utils for save 
def save_solution_to_file(file_path: str,
                 model_type: str,
                 model_value: float,
                 header: str, 
                 states: np.ndarray,) -> None:
    """
    Save the conditions to a csv file.

    The fist line presents the model propertie, folowed by the header and the data.

    MODEL_<model specific>=<model value>
    MODEL_IDT_T_400 -> IDT model, with T flag as 400 K
    MODEL_LFS -> LFS model
    """
    data_to_save = f"MODEL_{model_type}={model_value}\n{header}\n"
    data_to_save += "\n".join([",".join([str(states[row, col]) for col in range(states.shape[1])]) for row in range(states.shape[0])])
    
    with open(file_path, 'w') as f:
        f.write(data_to_save)


def create_solution_file_name(model_type: str, temp: float, press: float, phi: float) -> str:
    """
    Creates the file name considering the initial state.

    <model_type>_<temp>_<press>_<phi>.csv

    temp -> Kelvin
    press -> Pa

    if no phi is provided, use 1.0 as default.
    """
    f_name = f"{model_type}_{temp:.1f}_{press/ct.one_atm:.2f}_"
    if phi:
        f_name += f"{phi:.1f}"
    else:
        f_name += "1.0"
    f_name += ".csv"
    return f_name


# --------------------------------------------------------------------------------------------------
# IDT implementation
class ConstantMassPressureODE:
    """Implement the 0D, constant mass, constant pressure reactor"""
    def __init__(self, gas: ct.Solution) -> None:
        self.gas: ct.Solution = gas
        self._pressure: float = gas.P

    def __call__(self, t: float, y: np.ndarray) -> np.ndarray:
        """return the set of EDO. See Turns to understand."""
        # Set the gas conditions
        if y[0] <= 0:
            raise ValueError(f"Negative value found for temperature.")
        self.gas.set_unnormalized_mass_fractions(y[1:])
        self.gas.TP = y[0], self._pressure

        # calculates all the values
        rho = self.gas.density
        wdot = self.gas.net_production_rates
        dTdt = - (np.dot(self.gas.partial_molar_enthalpies, wdot) / (rho * self.gas.cp))
        dYdt = wdot * self.gas.molecular_weights / rho
        return np.hstack((dTdt, dYdt))


class IDTReactor:
    """
    Class implementation of an 0D reactor to obtain the IDT value.
    """
    def __init__(self, rtol: float = 10**-6, atol: float = 10**-9, n_iter: int = 10000) -> None:
        self.rtol = rtol
        self.atol = atol
        self.n_iter = n_iter
    
    def _get_solver(self, gas: ct.Solution, inital_condition: np.ndarray) -> integrate.ode:
        """
        Set the solver to run the IDT cases.
        """
        ODE = ConstantMassPressureODE(gas)
        solver = integrate.ode(ODE)
        solver.set_integrator('lsoda', method='bdf', rtol=self.rtol, atol=self.atol, with_jacobian=False, nsteps=self.n_iter)    
        solver.set_initial_value(inital_condition, 0.0)
        return solver

    def _get_header(self, gas: ct.Solution) -> str:
        return "time(s),T,P," + ",".join(["X_" + spc.name for spc in gas.species()])
        
    def get_idt(self, gas: ct.Solution, 
                    max_time: float = 5.0, 
                    dt: float = 10**-7, 
                    flag: str = 'T', 
                    temp_rise: float = 400.0) -> float:
        """
        Find the idt time.

        This returns a float. If an error is raised by solver problems, it returns a -1.0,
        if no IDT is found in the time window, returns -2.0.
        """
        # TODO: Add the species flag. For a peak, run up to the max time.

        # Prepare the initial conditions
        initial_condition = np.array([gas.T, *gas.Y])
        temperature_flag = temp_rise + initial_condition[0]

        # Set the solver
        solver = self._get_solver(gas, initial_condition)
        
        # solve
        try:
            while solver.successful() and solver.t <= max_time:
                if solver.y[0] >= temperature_flag:
                    return solver.t             
                solver.integrate(solver.t + dt)

        # catch any temperature problem
        except:
            return -1.0

        # if we do not find a IDT in the max_time
        return -2.0

    def get_norm_states(self, gas: ct.Solution, idt: float,
                            norm_dt: float = 0.01, 
                            max_norm_time: float = 2.0) -> Tuple[np.ndarray, str]:
        """
        Solve the idt reactor at every norm dt and return the real time conditions.
        
        Returns a np.ndarray containing the values and a str containig the header:
        time(s), T, P, Y_<spc names>
        """
        initial_condition = np.array([gas.T, *gas.Y])
        const_pressure = copy(gas.P)
        n_points = int(max_norm_time / norm_dt + 1)
        out_solution = np.zeros([n_points, 3 + gas.n_species])
        out_solution[0, 1], out_solution[0, 2]  = gas.T, gas.P
        out_solution[0, 3:] = gas.X

        # Set the solver
        solver = self._get_solver(gas, initial_condition)

        # set control parameters
        dt = norm_dt * idt
        max_time = max_norm_time * idt
        cur_point = 0

        try:
            while solver.successful():
                cur_point += 1
                solver.integrate(solver.t + dt)

                # for the output to be in molar fraction
                gas.TPY = solver.y[0], const_pressure, solver.y[1:]
                out_solution[cur_point, 0] = solver.t
                out_solution[cur_point, 1] = gas.T
                out_solution[cur_point, 2] = gas.P
                out_solution[cur_point, 3:] = gas.X
                
                if cur_point == n_points - 1:
                    break

            return out_solution, self._get_header(gas)
        except:
            raise Exception("Failed to solve the ODE. Try a different set of tolerances.")
    
    def solve(self, gas: ct.Solution, 
                    path_to_save: str = "",
                    phi_for_file_name: float = None,
                    norm_dt: float = 0.01, 
                    max_norm_time: float = 2.0,
                    max_time_for_idt: float = 5.0,
                    dt_for_idt: float = 10**-7,
                    idt_flag: str = 'T', 
                    idt_temp_rise: float = 400.0) -> Tuple[float, np.ndarray, str]:
        """
        Solve the reactor, returning a IDT value, a np.ndarray with all conditions and the corresponding header.
        If a directory is passed as input, save the conditions to the file.

        The condition in the gas solution passed to this method is consdered the initial condition.
        """
        init_TPY = copy(gas.TPY)
        idt = self.get_idt(gas, max_time=max_time_for_idt, dt=dt_for_idt, flag=idt_flag, temp_rise=idt_temp_rise)

        if idt <= 0.0:
            if idt == -2.0:
                raise Exception(f"It was not possble to obtain IDT. No IDT found in the time window.")
            raise Exception(f"It was not possble to obtain IDT. Solver problems found.")

        gas.TPY = init_TPY
        states, header = self.get_norm_states(gas, idt, norm_dt=norm_dt, max_norm_time=max_norm_time)

        # check for save flag:
        if path_to_save != "":
            f_name = create_solution_file_name("IDT", init_TPY[0], init_TPY[1], phi_for_file_name)
            save_solution_to_file(path_to_save + f_name, f"IDT_{idt_flag}_{idt_temp_rise:.2f}", idt, header, states)

        return idt, states, header
    

# --------------------------------------------------------------------------------------------------
# LFS implementation
class LFSReactor:
    """
    Class implementation of an 1D reactor to obtain the LFS value.
    """
    def __init__(self, width: float = 0.014, 
                        ratio: float=3, 
                        slope: float=0.1, 
                        curve: float=0.1,
                        max_time_step_count: int = 5000,
                        loglevel: int = 0) -> None:
        self.width = width
        self.ratio = ratio
        self.slope = slope
        self.curve = curve
        self.max_time_step_count = max_time_step_count
        self.loglevel = loglevel

    def _get_header(self, gas: ct.Solution) -> str:
        return "grid(s),T,P," + ",".join(["X_" + spc.name for spc in gas.species()])

    def _get_states(self, flame: ct.FreeFlame) -> np.ndarray:
        """
        Return the states of the current flame.

        grid(m), T, P, X_<species>
        """
        out_data = np.zeros([len(flame.T),len(flame.X) + 3])
        out_data[:,0] = flame.grid
        out_data[:,1] = flame.T
        out_data[:,2] = flame.P
        for i in range(len(flame.X)):
            out_data[:,3 + i] = flame.X[i]
        return out_data
    
    def solve(self, gas: ct.Solution, 
                    path_to_save: str = "",
                    phi_for_file_name: float = None) -> Tuple[float, np.ndarray, str]:
        """
        Solve the reactor, returning a IDT value, a np.ndarray with all conditions and the corresponding header.
        If a directory is passed as input, save the conditions to the file.

        The condition in the gas solution passed to this method is consdered the initial condition.
        """
        init_TPY = copy(gas.TPY)

        # Create the flame object
        flame = ct.FreeFlame(gas, width=self.width)
        # flame.transport_model = 'Mix'

        # Define tolerances for the solver
        flame.set_refine_criteria(ratio=self.ratio, slope=self.slope, curve=self.curve)
        flame.max_time_step_count = self.max_time_step_count

        # Define logging level
        flame.solve(loglevel=self.loglevel, auto=True)
        Su0 = flame.velocity[0]
        states = self._get_states(flame)
        header = self._get_header(gas)

        # check for save flag:
        if path_to_save != "":
            f_name = create_solution_file_name("LFS", init_TPY[0], init_TPY[1], phi_for_file_name)
            save_solution_to_file(path_to_save + f_name, f"LFS_m_s", Su0, header, states)

        return Su0, states, header






def main():
    gas = ct.Solution("gri30.xml")
    fuel = "CH4"
    oxi = "O2:1.0,N2:3.76"
    gas.set_equivalence_ratio(1.0, fuel, oxi)
    gas.TP = 1500.0, ct.one_atm
    dir_to_save = "C:\\Users\\1511 IRON\\Desktop\\Pós\\Doutorado\\PapersToBe\\Virtual_N_HEP\\red_teste4\\"

    # r = LFSReactor()
    # Su0 = r.solve(gas, path_to_save=dir_to_save)

    r = IDTReactor()
    # idt = r.get_idt(gas)
    # print(idt)
    idt, states, header = r.solve(gas, path_to_save=dir_to_save)


if __name__ == "__main__":
    main()