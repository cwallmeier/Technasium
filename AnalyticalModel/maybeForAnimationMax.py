"""
======================================================
@Project -> File: Technasium -> maybeForAnimationMax
@IDE            : PyCharm
@Author         : cwallmeier
@Date           : 29 Aug 25 22:20
======================================================
"""
import numpy as np
import pandas as pd
from doublet import geoth_doub
import matplotlib.pyplot as plt


# well spacing:
ws = 1000
# make a doublet
dub = geoth_doub(w_space=ws)

# A look at how the coldfront moves. Use the velocity function at a line of points from inj to prd
# give 10m buffer because the solution is a bit stupid near the wells.

x = pd.Series(np.linspace(10, ws-10, 50))
# speed of the front
v_cold = x.apply(dub.v_coldfront)

plt.plot(x, v_cold)
plt.ylim(0, max(v_cold) * 1.05)
plt.show()

# Tbh, this plot shows a very angular u-curve.. I think going wiht a constant average velocity for plotting purposes
