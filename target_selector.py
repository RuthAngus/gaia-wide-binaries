import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rotation import get_bv_and_age

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)

# load kepler-tgas catalogues.
star1_kic = pd.read_csv("star1_periods.csv")
star2_kic = pd.read_csv("star2_periods.csv")
star1_kic = get_bv_and_age(star1_kic)
star2_kic = get_bv_and_age(star2_kic)

# Plot 1: RA and DEC.
N = 150
plt.clf()
y1, x1 = star1_kic.dec_x.values[:N], star1_kic.ra_x.values[:N]
y2, x2 = star2_kic.dec_x.values[:N], star2_kic.ra_x.values[:N]
yerr1 = star1_kic.dec_error_x.values[:N]
yerr2 = star2_kic.dec_error_x.values[:N]
xerr1 = star1_kic.ra_error_x.values[:N]
xerr2 = star2_kic.ra_error_x.values[:N]
for i, ra in enumerate(x1):
    plt.plot([x1[i], x2[i]], [y1[i], y2[i]], color="k", alpha=.5)
plt.errorbar(x1, y1, yerr=yerr1, xerr=xerr1, fmt="k.", capsize=0, ecolor=".7")
plt.errorbar(x2, y2, yerr=yerr2, xerr=xerr2, fmt="k.", capsize=0, ecolor=".7")
plt.xlabel("$\mathrm{Right~Ascension~(degrees)}$")
plt.ylabel("$\mathrm{Declination~(degrees)}$")
plt.savefig("MDM_proposal/ra_vs_dec")

# # Plot 2: Rotation period vs colour?
# plt.clf()
# for i, age in enumerate(star1_kic.gyro_age.values):
#     if star1_kic.prot.values[i] > 0 and star2_kic.prot.values[i] > 0 \
#             and star1_kic.gyro_age.values[i] > 0 and \
#             star2_kic.gyro_age.values[i] > 0:
#             x = [star1_kic.B_V[i], star2_kic.B_V[i]]
#             y = [star1_kic.prot[i], star2_kic.prot[i]]
#             yerr = [star1_kic.prot_err[i], star2_kic.prot_err[i]]

#             plt.clf()
#             ys1 = period_model_b(star1_kic.gyro_age.values[i], xs)
#             ys2 = period_model_b(star2_kic.gyro_age.values[i], xs)
#             lab = star1_kic.gyro_age.values[i]
#             plt.plot(xs, ys1, "--", color="r",
#                     label="Age = {0:.3} Gyr".format(lab))
#             lab = star2_kic.gyro_age.values[i]
#             plt.plot(xs, ys2, "--", color="b",
#                     label="Age = {0:.3} Gyr".format(lab))
#             plt.errorbar(x[0], y[0], yerr=yerr[0], fmt="k.", ecolor=".7",
#                          capsize=0, markersize=12)
#             plt.errorbar(x[1], y[1], yerr=yerr[1], fmt="k.", ecolor=".7",
#                          capsize=0, markersize=12)

#             plt.ylabel("$\mathrm{Rotation~period~(days)}$")
#             plt.xlabel("$B-V$")
#             plt.legend()
#             plt.savefig("{0}_{1}_gyrochrones"
#                         .format(star1_kic.source_id.values[i],
#                                 star2_kic.source_id.values[i]))
