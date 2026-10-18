from itertools import cycle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import TABLEAU_COLORS, same_color
import typing
import scipy.stats as spt

"""
Calclulate and display a constraint diagram for the XTRA 2027 challenge ruleset
"""

class aircraft:
    """Holder class for aircraft containing all properties and from which performance is calculated"""

    def __init__(
        self,
        CL_max: float = 1.2,
        C_D0: float = 0.075,
        roll_resistance: float = 0.1, 
        TO_mass: float = 5, #kg
        AR: float = 10,
        oswald_e: float = 0.8,
        resolution: int = 20
    ) -> None:

        #params are separated into aircraft dependent and performance dependent
        #aircraft params
        self.CL_max = CL_max  # dimensionless
        self.C_D0 = C_D0
        self.roll_resistance = roll_resistance  # friction for take off
        self.g = 9.81  # gravity
        self.TO_mass = TO_mass  # kg
        self.weight = self.TO_mass * self.g  # N
        self.AR = AR
        self.oswald_e = oswald_e

        #find induced drag constant k from known aproximate properties
        self.K = 1 / (np.pi * self.oswald_e * self.AR)

        # plotting params
        self.resolution = resolution  # resolution for x scale on axis. W/S number of points
        self.ws = np.linspace(
            5, 90, self.resolution
        ) # construct array of take wing loading values to plot against

        #intialise an empty dictionary container for constraint results
        self.constraints = {}

    def _dynamic_p(self, v: float):
        """Return the dynamic pressure at sea level for a given velocity"""
        return 0.5 * 1.225 * v**2

    def ws_stall(self, v_stall: float):
        "Get the wing loading based on stall speed requirements"
        return self._dynamic_p(v_stall) * self.CL_max

    def tw_ground_roll(self, S_g: float, CD_lof: float, CL_lof: float, v_stall: float) -> np.ndarray:
        v_lof = 1.2 * v_stall
        q_lof = self._dynamic_p(0.7 * v_lof)

        tw = (
            (v_lof**2 / (2 * self.g * S_g))
            + (q_lof * CD_lof) / self.ws
            + self.roll_resistance * (1 - (q_lof * CL_lof) / self.ws)
        )

        return tw

    def tw_max_cruise_speed(self, v_cruise_max: float):
        """Calculate the T/W required for a given max cruise speed"""
        q_cvmax = self._dynamic_p(v_cruise_max)

        tw = ((q_cvmax * self.C_D0) / self.ws) + (self.K/q_cvmax) * self.ws

        return tw

    def _turn_load_factor_squared(self, v_turn: float, R_turn: float):
            """Helper to get the square of the load factor in a turn given the velocity and radius"""

            n_squared = (((v_turn ** 2) / (self.g * R_turn)) ** 2) + 1

            return n_squared

    def tw_turn_radius(self, v_turn, R_turn):
        """Calculate the T/W for a given turn radius at given turn speed"""
        q_turn = self._dynamic_p(v_turn)
        n_squared = self._turn_load_factor_squared(v_turn = v_turn, R_turn = R_turn)
        tw = ((q_turn * self.C_D0) / self.ws) + (self.K/q_turn) * self.ws * n_squared

        return tw

    def turn_stall_limit(self, v_turn, R_turn):
        """The maximum W/S for a sustained turn at given speed"""
        n = np.sqrt(self._turn_load_factor_squared(v_turn = v_turn, R_turn = R_turn))
        ws = self._dynamic_p(v = v_turn) * self.CL_max / n

        return ws

    def build_all_constraints(
            self,
            v_stall: float,
            S_g: float,
            v_cruise_max: float,
            R_turn: float,
            v_turn: float | None = None,
            CD_lof: float | None = None,
            CL_lof: float | None = None
    ):
        """The main builder function that computes all constraints and gives takes all target performace metrics"""

        #Set CD_lof to basic value based on cruise config with losses for low reynods if no entry
        if not CD_lof:
            CD_lof = self.C_D0*1.1

        #set CL_lof to crusie value if not provided. Account for low reynolds with scaling factor
        if not CL_lof:
            CL_lof = self.CL_max * 0.95

        #If no turn radius provided then set it the max speed
        if not v_turn:
            v_turn = v_cruise_max

        #Get various constraint tw values from entered list of target performance metrics
        #Key is the constraint name and the data type, must be array or value to be plotted correctley
        self.constraints.update(
            {
                ("ground_roll", "array"): self.tw_ground_roll(S_g = S_g, CD_lof = CD_lof, CL_lof = CL_lof, v_stall = v_stall),
                ("Stall_speed", "value"): self.ws_stall(v_stall = v_stall),
                ("max_cruise_speed", "array"): self.tw_max_cruise_speed(v_cruise_max = v_cruise_max),
                ("turn_radius", "array"): self.tw_turn_radius(v_turn = v_turn, R_turn = R_turn),
                ("turn_stall", "value"): self.turn_stall_limit(v_turn = v_turn, R_turn = R_turn)
            }
        )

    #builds graph with formatting
    def construct_graph(self):
        #set style of graphs
        plt.style.use('tableau-colorblind10')

        #create graph objects
        fig, ax = plt.subplots()

        #graph limits
        ax.set_ylim(0.1, 0.8) #hard limit on t/w for a reasonable range
        ax.set_xlim(self.ws[0], self.ws[-1])  #limit x axis to the computed range of w/s values 

        #Configure text displays on the graph
        ax.set_xlabel(r"$W/S$")
        ax.set_ylabel(r"$T/W$")
        ax.set_title("Sizing Constraint Diagram")

        #set colors
        colors = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

        #enable grid
        ax.grid(True, "both")

        return (fig, ax, colors)
    
    def plot_constraint_diagram(self):
        # Final output for the class, displays the constarint diagrams
        fig, ax, colors = self.construct_graph()

        #iterate over dictionary of constraints and plot
        for key, value in self.constraints.items():
            color = next(colors)
            #Choose whether to pot claculated vales as an array aginst w/s or as a vertical w/s line
            if key[1] == 'array':
                ax.plot(self.ws, value, label = key[0], color=color)
            elif key[1] == 'value':
                ax.vlines(value, 0.1, 0.5, label = key[0], color=color)
            else:
                raise ValueError(f"The passed data type is not supported: {value}")

        legend = ax.legend(loc = 'upper center')
        plt.show()


#Run program
def main():
    #run analysis for plane 1
    plane1 = aircraft()
    plane1.build_all_constraints(v_stall = 10, S_g = 40, v_cruise_max = 22, R_turn = 100)
    plane1.plot_constraint_diagram()


if __name__ == "__main__":
    main()