import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

DATA_DIR = "."
RESULTS_DIR = "rotation"
LC_DIR = "/Users/ruthangus/.kplr/data/lightcurves"


def teff2bv(teff, logg, feh):
    # best fit parameters
    t = [-813.3175, 684.4585, -189.923, 17.40875]
    f = [1.2136, 0.0209]
    d1 = -0.294
    g1 = -1.166
    e1 = 0.3125
    return t[0] + t[1] * np.log10(teff) + t[2] * (np.log10(teff))**2 + \
        t[3] * (np.log10(teff))**3 + f[0] * feh + f[1] * feh**2 + d1 * feh * \
        np.log10(teff) + g1 * logg + e1 * logg * np.log10(teff)


def age_model_mh(p, bv):
    a, n, b, c = .407, .566, .325, .495  # MH
    return (p / (a * (bv - c)**b))**(1./n) / 1000


def age_model_b(p, bv):
    a, b, c, n = .7725, .601, .4, .5189
    return (p / (a * (bv - c)**b))**(1./n) / 1000


def period_model_mh(age, bv):
    a, n, b, c = .407, .566, .325, .495  # MH
    age *= 1000
    return age**n * a * (bv - c)**b


def period_model_b(age, bv):
    a, b, c, n = .7725, .601, .4, .5189
    age *= 1000
    return age**n * a * (bv - c)**b


def search_db(id, df_name, DATA_DIR):
    prot, prot_err, ref = 0, 0, 0
    d = pd.read_csv(os.path.join(DATA_DIR, df_name))
    m = np.array(d["KIC"]) == int(id)
    if len(np.array(d["KIC"])[m]):
        prot = float(np.array(d["period"])[m])
        prot_err = float(np.array(d["period_err"])[m])
        ref = df_name
        print("Rotation period from {0}: \
                {1:.2} +/- {2:.2} Days".format(ref, prot, prot_err))
    return prot, prot_err, ref


def search_tables(id):
    periods1, period_errs1, refs1 = None, None, None
    periods1, period_errs1, refs1 = search_db(id, "vansaders.txt", ".")
    if not periods1:
        periods1, period_errs1, refs1 = \
            search_db(id, "Table_1_Periodic.txt", ".")
    if not periods1:
        periods1, period_errs1, refs1 = \
            search_db(id, "chaplin_garcia.csv", ".")
    if not periods1:
        Rdata = pd.read_csv(os.path.join(DATA_DIR,
                                         "Table_2_Non_Periodic.txt"))
        m = np.array(Rdata["KID"]) == int(id)
        if len(Rdata["KID"][m]):
            print("{0} is a non-rotating star".format(id))
            periods1, period_errs1, refs1 = 0, 0, \
                "Table_2_Non_Periodic.txt"
    return periods1, period_errs1, refs1


def get_periods(df):
    kids, periods, period_errs = [np.zeros(len(df)) for i in range(3)]
    refs = []
    for i, id in enumerate(df.kepid.values):
        print(id)
        periods[i], period_errs[i], ref = search_tables(id)
        refs.append(ref)
    return kids, periods, period_errs, refs


if __name__ == "__main__":

    # # load the kepler stars
    # star1_kic = pd.read_csv("star1_kic.csv")
    # star2_kic = pd.read_csv("star2_kic.csv")

    # kids1, periods1, period_errs1, refs1 = get_periods(star1_kic)
    # kids2, periods2, period_errs2, refs2 = get_periods(star2_kic)

    # # save periods and refs
    # print(periods1)
    # star1_kic["prot"] = periods1
    # star1_kic["prot_err"] = period_errs1
    # star1_kic["prot_ref"] = refs1
    # star1_kic.to_csv("star1_periods.csv")
    # star2_kic["prot"] = periods2
    # star2_kic["prot_err"] = period_errs2
    # star2_kic["prot_ref"] = refs2
    # star2_kic.to_csv("star2_periods.csv")

    star1_kic = pd.read_csv("star1_ages.csv")
    star2_kic = pd.read_csv("star2_ages.csv")

    xs = np.arange(.4, 1.5, .01)
    for i, age in enumerate(star1_kic.gyro_age.values):
        if star1_kic.prot.values[i] > 0 and star2_kic.prot.values[i] > 0:
            x = [star1_kic.B_V[i], star2_kic.B_V[i]]
            y = [star1_kic.prot[i], star2_kic.prot[i]]
            yerr = [star1_kic.prot_err[i], star2_kic.prot_err[i]]

            plt.clf()
            plt.errorbar(x[0], y[0], yerr=yerr[0], fmt="r.")
            plt.errorbar(x[1], y[1], yerr=yerr[1], fmt="b.")
            plt.plot(x, y, color=".7")
            ys1 = period_model_b(star1_kic.gyro_age.values[i], xs)
            ys2 = period_model_b(star2_kic.gyro_age.values[i], xs)
            plt.plot(xs, ys1, "--", color="r",
                     label="{0:.3} Gyr, log(g) = "
                     "{1:.3}, {2}, {3}".format(star1_kic.gyro_age[i],
                                               star1_kic.logg[i],
                                               star1_kic.prot_ref[i],
                                               star1_kic.kepid[i]))
            plt.plot(xs, ys2, "--", color="b",
                     label="{0:.3} Gyr, log(g) = "
                     "{1:.3}, {2}, {3}".format(star2_kic.gyro_age[i],
                                               star2_kic.logg[i],
                                               star2_kic.prot_ref[i],
                                               star2_kic.kepid[i]))
            plt.ylabel("period")
            plt.xlabel("B-V")
            plt.legend()
            plt.savefig("{0}".format(str(i).zfill(2)))
