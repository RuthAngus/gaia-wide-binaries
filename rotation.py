import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import h5py
import glob

plotpar = {'axes.labelsize': 18,
           'font.size': 10,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)


DATA_DIR = "data"
RESULTS_DIR = "rotation"
LC_DIR = "/Users/ruthangus/.kplr/data/lightcurves"


def teff2bv(teff, logg, feh):
    """
    Convert teff to B-V colour.
    """
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
    """
    Mamajek & Hillenbrand (2008) age model.
    """
    a, n, b, c = .407, .566, .325, .495  # MH
    return (p / (a * (bv - c)**b))**(1./n) / 1000


def age_model_b(p, bv):
    """
    Barnes (2003) age model.
    """
    a, b, c, n = .7725, .601, .4, .5189
    return (p / (a * (bv - c)**b))**(1./n) / 1000


def period_model_mh(age, bv):
    """
    Mamajek & Hillenbrand (2008) period model.
    """
    a, n, b, c = .407, .566, .325, .495  # MH
    age *= 1000
    return age**n * a * (bv - c)**b


def period_model_b(age, bv):
    """
    Barnes (2003) period model.
    """
    a, b, c, n = .7725, .601, .4, .5189
    age *= 1000
    return age**n * a * (bv - c)**b


def search_db(id, df_name, DATA_DIR):
    """
    Search a database (text file) a kic id and find the rotation period.
    id: Kepid
    df_name: the path to the text file. E.g. 'Table_1_Periodic.txt' for
    McQuillan, 'data/chaplin-garcia.csv' for garcia or 'data/vansaders.txt'
    for van saders.
    """
    prot, prot_err, ref = 0., 0., 0.
    d = pd.read_csv(os.path.join(DATA_DIR, df_name))
    m = d.KIC.values == int(id)
    if len(d.KIC.values[m]):
        prot = float(d.period.values[m])
        prot_err = float(d.period_err.values[m])
        ref = df_name
        print("Rotation period from {0}: \
                {1:.2} +/- {2:.2} Days".format(ref, prot, prot_err))
    return prot, prot_err, ref


def search_tables(id, DATA_DIR):
    """
    Search several tables for rotation periods.
    id: kepid
    DATA_DIR: the directory of the rotation period catalogue file.
    """
    periods1, period_errs1, refs1 = search_db(id, "Table_1_Periodic.txt",
                                              DATA_DIR)
    if not periods1:
        periods1, period_errs1, refs1 = \
            search_db(id, "data/chaplin_garcia.csv", DATA_DIR)
    if not periods1:
        periods1, period_errs1, refs1 = \
            search_db(id, "data/vansaders.txt", DATA_DIR)
    return periods1, period_errs1, refs1


def get_periods(df):
    """
    Find a period and reference for that period from the literature for all
    stars in dataframe, df.
    """
    kids, periods, period_errs = [np.zeros(len(df)) for i in range(3)]
    refs = []
    for i, id in enumerate(df.kepid.values):
        print(id)
        periods[i], period_errs[i], ref = search_tables(id)
        refs.append(ref)
    return kids, periods, period_errs, refs


def get_bv_and_age(df):
    """
    Calculate the B-V colour and age of the stars in dataframe df.
    """
    # estimate B-V
    df["B_V"] = teff2bv(df["teff"], df["logg"], df["feh"])

    ages = np.zeros_like(df.B_V.values)
    m = df.B_V.values > 0.4
    ages[m] = age_model_b(df["prot"][m], df["B_V"][m])
    df["gyro_age"] = ages

    fname = "data/chaplin_garcia.csv"
    dat = pd.read_csv(fname)
    astero_ages = np.zeros_like(df.prot.values)
    astero_age_errps = np.zeros_like(df.prot.values)
    astero_age_errms = np.zeros_like(df.prot.values)
    for i, _ in enumerate(df.prot_ref):
        if df.prot_ref.values[i] == fname:
            m = dat["KIC"] == df.kepid[i]
            astero_ages[i] == dat["age"]
            astero_age_errps[i] == dat["age_errp"]
            astero_age_errms[i] == dat["age_errm"]
    df["astero_age"] = astero_ages
    df["astero_age_errp"] = astero_age_errps
    df["astero_age_errm"] = astero_age_errps
    return df


def find_wbstar_prots():
    # load the kepler stars
    star1_kic = pd.read_csv("star1_kic.csv")
    star2_kic = pd.read_csv("star2_kic.csv")

    kids1, periods1, period_errs1, refs1 = get_periods(star1_kic)
    kids2, periods2, period_errs2, refs2 = get_periods(star2_kic)

    # save periods and refs
    print(periods1)
    star1_kic["prot"] = periods1
    star1_kic["prot_err"] = period_errs1
    star1_kic["prot_ref"] = refs1
    star1_kic.to_csv("star1_periods.csv")
    star2_kic["prot"] = periods2
    star2_kic["prot_err"] = period_errs2
    star2_kic["prot_ref"] = refs2
    star2_kic.to_csv("star2_periods.csv")


def plot_binary_gyrochrones():

    star1_kic = pd.read_csv("star1_periods.csv")
    star2_kic = pd.read_csv("star2_periods.csv")

    star1_kic = get_bv_and_age(star1_kic)
    star2_kic = get_bv_and_age(star2_kic)

    xs = np.arange(.4, 1.5, .01)
    for i, age in enumerate(star1_kic.gyro_age.values):

        if star1_kic.prot.values[i] > 0 and star2_kic.prot.values[i] > 0 \
                and star1_kic.gyro_age.values[i] > 0 and \
                star2_kic.gyro_age.values[i] > 0:
            x = [star1_kic.B_V[i], star2_kic.B_V[i]]
            y = [star1_kic.prot[i], star2_kic.prot[i]]
            yerr = [star1_kic.prot_err[i], star2_kic.prot_err[i]]

            plt.clf()
            # plt.plot(x, y, color=".7")
            ys1 = period_model_b(star1_kic.gyro_age.values[i], xs)
            ys2 = period_model_b(star2_kic.gyro_age.values[i], xs)
            lab = star1_kic.gyro_age.values[i]
            plt.plot(xs, ys1, "--", color="r",
                    label="Age = {0:.3} Gyr".format(lab))
            lab = star2_kic.gyro_age.values[i]
            plt.plot(xs, ys2, "--", color="b",
                    label="Age = {0:.3} Gyr".format(lab))
            plt.errorbar(x[0], y[0], yerr=yerr[0], fmt="k.", ecolor=".7",
                         capsize=0, markersize=12)
            plt.errorbar(x[1], y[1], yerr=yerr[1], fmt="k.", ecolor=".7",
                         capsize=0, markersize=12)

            plt.ylabel("$\mathrm{Rotation~period~(days)}$")
            plt.xlabel("$B-V$")
            plt.legend()
            plt.savefig("{0}_{1}_gyrochrones"
                        .format(star1_kic.source_id.values[i],
                                star2_kic.source_id.values[i]))


def posteriors():

    star1_kic = pd.read_csv("star1_periods.csv")
    star2_kic = pd.read_csv("star2_periods.csv")
    print(star1_kic.keys())
    print(star1_kic.ra_x.values[0], star1_kic.dec_x.values[0])

    star1_kic = get_bv_and_age(star1_kic)
    star2_kic = get_bv_and_age(star2_kic)

    fnames = glob.glob("samples/*h5")

    med_lnage0, med_lnage1, std_lnage0, std_lnage1 = [], [], [], []
    g_age0, g_age1 = [], []
    for fname in fnames:
        si1 = fname[8:27]
        si2 = fname[28:-3]
        m1 = star1_kic.source_id.values == int(si1)
        m2 = star2_kic.source_id.values == int(si2)
        print(star1_kic.source_id.values[m1][0])
        print(star2_kic.source_id.values[m2][0])

        gyro_age1 = star1_kic.gyro_age.values[m1][0]
        gyro_age2 = star2_kic.gyro_age.values[m2][0]

        # get age posterior
        samples = pd.read_hdf(fname, "df")
        med_lnage0.append(np.median(samples.age_0_0.values))
        med_lnage1.append(np.median(samples.age_0_1.values))

        # plt.clf()
        # # plt.hist(samples.age_0_0.values, color="r", alpha=.5)
        # plt.hist(samples.age_0_1.values, histtype="stepfilled", color="w",
        #          label="{0:.2} Gyr".format(age1))
        # plt.axvline(np.log10(gyro_age1 * 1e9), color="r", ls="--",
        #             label="{0:.2} Gyr".format(gyro_age1))
        # plt.axvline(np.log10(gyro_age2 * 1e9), color="b", ls="--",
        #             label="{0:.2} Gyr".format(gyro_age2))
        # plt.legend(loc="best")
        # plt.xlabel("$\ln(Age)$")
        # plt.savefig("{0}_{1}_hist".format(si1, si2))

        std_lnage0.append(np.std(samples.age_0_0.values))
        std_lnage1.append(np.std(samples.age_0_1.values))
        g_age0.append(gyro_age1)
        g_age1.append(gyro_age2)

        # plt.clf()
        # x = range(13)
        # plt.errorbar(x, np.median(samples.age_0_1.values), yerr=yerr,
        #             fmt="ko")
        # # plt.plot(x, np.median(samples.age_0_1.values), "ko", alpha=.7)
        # plt.plot(x, np.log10(gyro_age1 * 1e9), "ro", alpha=.7)
        # plt.plot(x, np.log10(gyro_age2 * 1e9), "bo", alpha=.7)
        # plt.xlabel("$\mathrm{Star~index}$")
        # plt.ylabel("$\ln(\mathrm{Age})$")
        # plt.xlim(1, 13)
        # plt.savefig("age_summary")

    plt.clf()
    x = range(1, 13)
    plt.errorbar(x, med_lnage0, yerr=std_lnage0, fmt="ko")
    # plt.plot(x, np.median(samples.age_0_1.values), "ko", alpha=.7)
    plt.plot(x, np.log10(np.array(g_age0) * 1e9), "ro", alpha=.7)
    plt.plot(x, np.log10(np.array(g_age1) * 1e9), "bo", alpha=.7)
    plt.xlabel("$\mathrm{Star~index}$")
    plt.ylabel("$\log_{10}(\mathrm{Age})$")
    plt.xlim(0, 13)
    plt.savefig("age_summary")
    print("1")

    x = np.arange(1, 13)
    # mm = np.ones_like(x)
    # mm[2], mm[5], mm[10] = 0, 0, 0
    # m = mm == 1
    m = x < 50

    orange = '#FF9933'
    lightblue = '#66CCCC'
    blue = '#0066CC'
    pink = '#FF33CC'
    turquoise = '#3399FF'
    lightgreen = '#99CC99'
    green = '#009933'
    maroon = '#CC0066'
    purple = '#9933FF'
    red = '#CC0000'
    lilac = '#CC99FF'

    print("2")
    med_lnage0 = np.array(med_lnage0)[m]
    med_lnage1 = np.array(med_lnage1)[m]
    std_lnage0 = np.array(std_lnage0)[m]
    std_lnage1 = np.array(std_lnage1)[m]
    g_age0 = np.array(g_age0)[m]
    g_age1 = np.array(g_age1)[m]
    inds = np.array(np.argsort(med_lnage0))

    med_lnage0 = med_lnage0[inds]
    med_lnage1 = med_lnage1[inds]
    std_lnage0 = std_lnage0[inds]
    std_lnage1 = std_lnage1[inds]
    g_age0 = g_age0[inds]
    g_age1 = g_age1[inds]


    plt.clf()
    # x = np.arange(1, len(med_lnage0)+1)
    print(len(x), len(med_lnage0), len(std_lnage0))
    plt.errorbar(x, med_lnage0, yerr=std_lnage0, fmt="ko", capsize=0,
                 ms=10, alpha=.7)
    # plt.plot(x, np.median(samples.age_0_1.values), "ko", alpha=.7)
    plt.plot(x, np.log10(g_age0 * 1e9), "o", alpha=.7, ms=10,
             mec="None", color=pink)
    plt.plot(x, np.log10(g_age1 * 1e9), "o", alpha=.7, ms=10,
             mec="None", color=blue)
    plt.xlabel("$\mathrm{Star~index}$")
    plt.ylabel("$\log_{10}(\mathrm{Age})$")
    # plt.xlim(0, 10)
    plt.savefig("age_summary3")

    med_age0 = 10**(np.array(med_lnage0))*1e-9
    med_age1 = 10**(np.array(med_lnage1))*1e-9
    std_age0 = 10**(np.array(std_lnage0)*np.array(med_age0))*1e-9
    std_age1 = 10**(np.array(std_lnage1))*1e-9

    plt.clf()
    x = range(1, 13)
    plt.errorbar(x, med_age0, yerr=std_age0, fmt="ko")
    # plt.plot(x, np.median(samples.age_0_1.values), "ko", alpha=.7)
    plt.plot(x, g_age0, "ro", alpha=.7)
    plt.plot(x, g_age1, "bo", alpha=.7)
    plt.xlabel("$\mathrm{Star~index}$")
    plt.ylabel("$\ln(\mathrm{Age})$")
    plt.xlim(0, 13)
    # plt.ylim(0, 6)
    plt.savefig("age_summary2")

if __name__ == "__main__":

    # plot_binary_gyrochrones()
    plot_binary_gyrochrones()
    assert 0
    posteriors()

    # find all tgas kepler rotation periods

    # load kepler-tgas
    kplr_tgas = pd.read_csv("kic_tgas.csv")

    # calculate B-Vs
    kplr_tgas["B_V"] = teff2bv(kplr_tgas.teff, kplr_tgas.logg, kplr_tgas.feh)

    # look up rotation periods
    prot, prot_err = [np.zeros_like(kplr_tgas.kepid.values) for i in range(2)]
    refs = []
    for i, id in enumerate(kplr_tgas.kepid.values):
        print(i, "of", len(kplr_tgas.kepid.values))
        p, p_err, ref = search_tables(id, DATA_DIR)
        prot[i], prot_err[i] = p, p_err
        refs.append(ref)
        with open("prots.txt", "a") as f:
            f.write("{0} {1} {2} {3} \n".format(id, p, p_err, ref))

    # Save periods and b-vs
    kplr_tgas["prot"] = prot
    kplr_tgas["prot_err"] = prot_err
    kplr_tgas["prot_ref"] = refs
    kplr_tgas.to_csv("kplr_tgas_periods.csv")
