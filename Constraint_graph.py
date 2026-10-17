import numpy as np
from matplotlib import pyplot as plt
import typing
import scipy.stats as spt

"""
Calclulate and display a constraint diagram for the XTRA 2027 challenge ruleset
"""

# Define hard aircraft params


class aircraft:
    """Holder class for aircraft containing all properties and from which performance is calculated"""

    def __init__(
        self,
        CL_max: float = 1.2,
        v_stall: float = 12,
        C_D0: float = 0.075,
        roll_resistance: float = 0.1,
        TO_mass: float = 5,
        resolution: int = 20,
    ) -> None:

        self.CL_max = CL_max  # dimensionless
        self.v_stall = v_stall  # m/s
        self.C_D0 = C_D0
        self.roll_resistance = roll_resistance  # friction for take off
        self.g = 9.81  # gravity
        self.TO_mass = TO_mass  # kg
        self.weight = self.TO_mass * self.g  # N

        # plotting params
        self.resolution = (
            resolution  # resolution for x scale on axis. W/S number of points
        )
        self.ws = np.linspace(
            10, 50, self.resolution
        ) # construct array of take wing loading values to plot against

    def _dynamic_p(self, v):
        """Return the dynamic pressure at sea level for a given velocity"""
        return 0.5 * 1.225 * v**2

    def ws_stall(self):
        "Get the wing loading based on stall speed requirements"
        return self._dynamic_p(self.v_stall) * self.CL_max

    def tw_ground_roll(self, S_g, CD_lof, CL_lof) -> np.ndarray:
        v_lof = 1.2 * self.v_stall
        q_lof = self._dynamic_p(0.7 * v_lof)

        tw = (
            (v_lof**2 / (2 * self.g * S_g))
            + (q_lof * CD_lof) / self.ws
            + self.roll_resistance * (1 - (q_lof * CL_lof) / self.ws)
        )

        return tw

    #This is to plot my diagram
    def plot_constraint_diagram(self):
        # Final output for the class, displays the constarint diagrams
        plot, axes = plt.Figure()
        axes.plot(...)
        return None


if __name__ == "__main__":
    plane1 = aircraft()
    tw_lift_off = plane1.tw_ground_roll(60, 0.05, 0.7)

    print(tw_lift_off)
    print(__name__ == "__main__")
