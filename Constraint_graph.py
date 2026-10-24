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
        CL_max: float,
        C_D0: float,
        roll_resistance: float,
        AR: float,
        oswald_e: float
    ) -> None:

        # params are separated into aircraft dependent and performance dependent
        # environemnt params
        self.rho: float = 1.225  # air density
        self.g: float = 9.81  # gravity

        # aircraft params
        self.CL_max = CL_max  # dimensionless
        self.C_D0 = C_D0
        self.roll_resistance = roll_resistance  # friction for take off
        self.AR = AR
        self.oswald_e = oswald_e

        # find induced drag constant k from known aproximate properties
        self.K = 1 / (np.pi * self.oswald_e * self.AR)

        # plotting params
        self.resolution = 45  # resolution for x scale on axis. W/S number of points

        self.ws = np.linspace(
            5, 150, self.resolution
        )  # construct array of take wing loading values to plot against

        # intialise an empty dictionary container for constraint results
        self.constraints = {}

    def _dynamic_p(self, v: float):
        """Return the dynamic pressure at sea level for a given velocity"""
        return 0.5 * self.rho * v**2

    def ws_stall(self, v_stall: float):
        "Get the wing loading based on stall speed requirements"
        return self._dynamic_p(v_stall) * self.CL_max

    def stall_velocity(self) -> np.ndarray:
        """Calculate stall velocity for a range of wing loadings."""
        return np.sqrt(2 * self.ws / (self.rho * self.CL_max))

    def tw_ground_roll(
        self, S_g: float, CD_lof: float, CL_lof: float, v_stall: float
    ) -> np.ndarray:
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

        tw = ((q_cvmax * self.C_D0) / self.ws) + (self.K / q_cvmax) * self.ws

        return tw

    def _turn_load_factor_squared(self, v_turn: float, R_turn: float):
        """Helper to get the square of the load factor in a turn given the velocity and radius"""
        n_squared = (((v_turn**2) / (self.g * R_turn)) ** 2) + 1

        return n_squared

    def tw_turn_radius(self, v_turn, R_turn):
        """Calculate the T/W for a given turn radius at given turn speed"""
        q_turn = self._dynamic_p(v_turn)
        n_squared = self._turn_load_factor_squared(v_turn=v_turn, R_turn=R_turn)

        tw = ((q_turn * self.C_D0) / self.ws) + (self.K / q_turn) * self.ws * n_squared

        return tw

    def turn_stall_limit(self, v_turn, R_turn):
        """The maximum W/S for a sustained turn at given speed"""
        n = np.sqrt(self._turn_load_factor_squared(v_turn=v_turn, R_turn=R_turn))
        ws = self._dynamic_p(v=v_turn) * self.CL_max / n

        return ws

    def _v_min_drag(self):
        "find the minimum drag velocity of the aircraft for rage of wing loadings"
        v_md = np.sqrt((2 * self.ws / self.rho) * np.sqrt(self.K / self.C_D0))

        return v_md

    def tw_climb_rate(self, climb_rate: float):
        """Get lisat of T/W values for a given climb rate"""
        v_md = self._v_min_drag()
        q_climb = self._dynamic_p(v_md)

        tw = (
            ((q_climb * self.C_D0) / self.ws)
            + ((self.K / q_climb) * self.ws)
            + (climb_rate / v_md)
        )

        return tw

    def build_all_constraints(
        self,
        v_stall: float,
        S_g: float,
        v_cruise_max: float,
        R_turn: float,
        climb_rate: float,
        v_turn: float | None = None,
        CD_lof: float | None = None,
        CL_lof: float | None = None
    ):
        """The main builder function that computes all constraints and gives takes all target performace metrics"""

        # Set CD_lof to basic value based on cruise config with losses for low reynods if no entry
        if not CD_lof:
            CD_lof = self.C_D0 * 1.1

        # set CL_lof to crusie value if not provided. Account for low reynolds with scaling factor
        if not CL_lof:
            CL_lof = self.CL_max * 0.95

        # If no turn radius provided then set it the max speed
        if not v_turn:
            v_turn = v_cruise_max

        # Store each constraint by name and the plot type used to display it.
        self.constraints.update(
            {
                ("ground_roll", "array"): self.tw_ground_roll(
                    S_g=S_g, CD_lof=CD_lof, CL_lof=CL_lof, v_stall=v_stall
                ),
                ("stall_speed", "secondary_array"): self.stall_velocity(),
                ("max_cruise_speed", "array"): self.tw_max_cruise_speed(
                    v_cruise_max=v_cruise_max
                ),
                ("turn_radius", "array"): self.tw_turn_radius(
                    v_turn=v_turn, R_turn=R_turn
                ),
                ("turn_stall", "value"): self.turn_stall_limit(
                    v_turn=v_turn, R_turn=R_turn
                ),
                ("climb_rate", "array"): self.tw_climb_rate(climb_rate=climb_rate),
            }
        )

    # builds graph with formatting
    def construct_graph(
        self,
        aircraft_params: dict[str, float | int],
        performance_targets: dict[str, float | int],
    ):
        plt.style.use("tableau-colorblind10")

        fig = plt.figure(figsize=(13, 7), layout="constrained")
        grid = fig.add_gridspec(2, 2, width_ratios=(3.2, 1.5))
        ax = fig.add_subplot(grid[:, 0])
        aircraft_info_ax = fig.add_subplot(grid[0, 1])
        performance_info_ax = fig.add_subplot(grid[1, 1])

        ax.set_xlim(self.ws[0], self.ws[-1])
        ax.set_ylim(0.1, 0.8)

        ax.set_xlabel(r"$W/S$ (N/m$^2$)", labelpad=10)
        ax.set_ylabel(r"$T/W$", labelpad=10)
        ax.set_title("Aircraft Performance Constraints", pad=14, weight="bold")

        colors = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

        ax.set_axisbelow(True)
        ax.minorticks_on()
        ax.grid(which="major", color="#b8bec5", linewidth=0.8, alpha=0.7)
        ax.grid(which="minor", color="#d9dde1", linewidth=0.5, linestyle=":", alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", which="major", direction="out", length=5)
        ax.tick_params(axis="both", which="minor", direction="out", length=3)

        aircraft_labels = {
            "CL_max": ("Maximum lift coefficient", ""),
            "C_D0": ("Zero-lift drag coefficient", ""),
            "roll_resistance": ("Takeoff rolling resistance", ""),
            "AR": ("Aspect ratio", ""),
            "oswald_e": ("Oswald efficiency factor", "")
        }
        performance_labels = {
            "v_stall": ("Stall speed", "m/s"),
            "S_g": ("Ground roll distance", "m"),
            "v_cruise_max": ("Maximum cruise speed", "m/s"),
            "climb_rate": ("Climb rate", "m/s"),
            "R_turn": ("Turn radius", "m"),
            "v_turn": ("Turn speed", "m/s"),
        }

        def format_parameters(
            parameters: dict[str, float | int],
            labels: dict[str, tuple[str, str]],
        ) -> str:
            lines = []
            for name, value in parameters.items():
                label, unit = labels.get(name, (name.replace("_", " ").title(), ""))
                formatted_value = f"{value:g}" if isinstance(value, (int, float)) else str(value)
                suffix = f" {unit}" if unit else ""
                lines.append(f"{label}: {formatted_value}{suffix}")
            return "\n".join(lines)

        for info_ax, title, parameters, labels in (
            (
                aircraft_info_ax,
                "Aircraft parameters",
                aircraft_params,
                aircraft_labels
            ),
            (
                performance_info_ax,
                "Assumed performance metrics",
                performance_targets,
                performance_labels,
            ),
        ):
            info_ax.set_axis_off()
            info_ax.text(
                0.03,
                0.96,
                title,
                transform=info_ax.transAxes,
                va="top",
                fontsize=11,
                fontweight="bold",
            )
            info_ax.text(
                0.03,
                0.84,
                format_parameters(parameters, labels),
                transform=info_ax.transAxes,
                va="top",
                fontsize=9,
                linespacing=1.6,
                bbox={
                    "boxstyle": "round,pad=0.7",
                    "facecolor": "white",
                    "edgecolor": "#b8bec5",
                    "linewidth": 0.8,
                },
            )

        return fig, ax, colors

    def plot_constraint_diagram(
        self,
        aircraft_params: dict[str, float | int],
        performance_targets: dict[str, float | int],
    ):
        # Final output for the class, displays the constarint diagrams
        fig, ax, colors = self.construct_graph(aircraft_params, performance_targets)
        stall_ax = ax.twinx()
        stall_ax.spines["top"].set_visible(False)
        stall_ax.set_ylabel("Stall velocity (m/s)")

        # iterate over dictionary of constraints and plot
        for key, value in self.constraints.items():
            color = next(colors)
            # Choose whether to pot claculated vales as an array aginst w/s or as a vertical w/s line
            if key[1] == "array":
                ax.plot(self.ws, value, label=key[0], color=color)

            elif key[1] == "secondary_array":
                stall_ax.plot(
                    self.ws,
                    value,
                    label=key[0],
                    color=color,
                    linestyle="dashed",
                    alpha = 0.8
                )

            elif key[1] == "value":
                ax.vlines(value, 0, 2, label=key[0], color=color, linestyles="dashed")

            elif key[1] == "point":
                raise NotImplementedError

            else:
                raise ValueError(f"The passed data type is not supported: {value}")

        handles, labels = ax.get_legend_handles_labels()
        stall_handles, stall_labels = stall_ax.get_legend_handles_labels()
        ax.legend(handles + stall_handles, labels + stall_labels, loc="upper right")
        plt.show()


# main function to instantiate classes and give aircraft params
def main():
    aircraft1_params = {
        "CL_max": 1.2,
        "C_D0": 0.09,
        "roll_resistance": 0.1,
        "AR": 12,
        "oswald_e": 0.75
    }
    aircraft1_performance_targets = {
        "v_stall": 10,
        "S_g": 40,
        "v_cruise_max": 20,
        "R_turn": 50,
        "climb_rate": 1.66,
        "v_turn": 20
    }

    # run analysis for plane 1
    plane1 = aircraft(**aircraft1_params)
    plane1.build_all_constraints(**aircraft1_performance_targets)
    plane1.plot_constraint_diagram(aircraft1_params, aircraft1_performance_targets)


# run program
if __name__ == "__main__":
    main()
