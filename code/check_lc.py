""" Check the 1DC LC of the surviving sources (not just programid=3)

Only keep sources that have,
- detections on three consecutive nights
- first detection fainter than peak
- last detection fainter than peak
"""

import numpy as np
import glob


