import numpy as np
import pyfits
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt


def xmatch(m1, m2, fname):
    """
    Match the stars to the Kepler and K2 catalogues.
    """
    # load full TGAS
    table = Table.read('stacked_tgas.fits')
    stacked_tgas_df = table.to_pandas()

    # get binaries in TGAS
    star1_df = stacked_tgas_df.iloc[m1]
    star2_df = stacked_tgas_df.iloc[m2]

    # load kepler-tgas xmatch
    kplr_tgas = pd.read_csv("kic_tgas_mod.csv")
    epic_tgas = pd.read_csv("epic_tgas_mod.csv")

    # find the indices of the star1s and 2s that are in kepler-tgas
    m_tgas = star1_df.source_id.isin(kplr_tgas.source_id).values * \
        star2_df.source_id.isin(kplr_tgas.source_id).values
    star1_tgas = star1_df.iloc[np.where(m_tgas)[0]]
    star2_tgas = star2_df.iloc[np.where(m_tgas)[0]]
    print(len(np.where(m_tgas)[0]), "pairs found in kepler")

    # find the indices of the star1s and 2s that are in epic-tgas
    m_tgas_epic = star1_df.source_id.isin(epic_tgas.source_id).values * \
        star2_df.source_id.isin(epic_tgas.source_id).values
    star1_tgas_epic = star1_df.iloc[np.where(m_tgas_epic)[0]]
    star2_tgas_epic = star2_df.iloc[np.where(m_tgas_epic)[0]]
    print(len(np.where(m_tgas_epic)[0]), "pairs found in K2")

    star1_kic = pd.merge(star1_tgas, kplr_tgas, on="source_id")
    star1_kic.to_csv("star1_kic_{0}.csv".format(fname))
    star2_kic = pd.merge(star2_tgas, kplr_tgas, on="source_id")
    star2_kic.to_csv("star2_kic_{0}.csv".format(fname))
    star1_epic = pd.merge(star1_tgas_epic, epic_tgas, on="source_id")
    star1_epic.to_csv("star1_epic_{0}.csv".format(fname))
    star2_epic = pd.merge(star2_tgas_epic, epic_tgas, on="source_id")
    star2_epic.to_csv("star2_epic_{0}.csv".format(fname))

    plt.clf()
    for i, _ in enumerate(star1_kic.ra_x.values):
        plt.plot([star1_kic.ra_x.values[i], star2_kic.ra_x.values[i]],
                 [star1_kic.dec_x.values[i], star2_kic.dec_x.values[i]])
    plt.plot(star1_kic.ra_x, star1_kic.dec_x, "r.")
    plt.plot(star2_kic.ra_x, star2_kic.dec_x, "k.")
    plt.savefig("kepler_ra_dec")

    plt.clf()
    for i, _ in enumerate(star1_epic.ra_x.values):
        plt.plot([star1_epic.ra_x.values[i], star2_epic.ra_x.values[i]],
                 [star1_epic.dec_x.values[i], star2_epic.dec_x.values[i]])
    plt.plot(star1_epic.ra_x, star1_epic.dec_x, "k.")
    plt.plot(star2_epic.ra_x, star2_epic.dec_x, "k.")
    # plt.xlim(230, 270)
    # plt.ylim(-30, -12)
    plt.savefig("k2_ra_dec")

if __name__ == "__main__":

    # load TGAS star indices
    hdulist = pyfits.open("pairindices_cp1.fits")
    t = hdulist[1].data
    m1 = t["star1"]
    m2 = t["star2"]
    xmatch(m1, m2, "1")
