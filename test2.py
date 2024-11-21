import os

import numpy as np
from astropy import units as u

from orbits import Orbit, interplanetary_transfer_dv, OrbitalElements
import two_body_util as util

width = os.get_terminal_size().columns

# Question 1:  ---------------------------
print("\n" + "-"*width + "\n")
print("Question 1: \n")
