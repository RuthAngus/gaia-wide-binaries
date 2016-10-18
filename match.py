import numpy as np
import pyfits
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt

# load TGAS star indices
hdulist = pyfits.open("pairindices_cp01.fits")
t = hdulist[1].data
m1 = t["star1"]
m2 = t["star2"]

# load full TGAS
hdulist_tgas = pyfits.open("stacked_tgas.fits")
table = Table.read('stacked_tgas.fits')
stacked_tgas_df = table.to_pandas()

# get binaries in TGAS
star1_df = stacked_tgas_df.iloc[m1]
star2_df = stacked_tgas_df.iloc[m2]

keep_track = np.vstack((star1_df.source_id, star2_df.source_id)).T

# load kepler-tgas xmatch
kplr_tgas = pd.read_csv("kic_tgas_mod.csv")
epic_tgas = pd.read_csv("epic_tgas_mod.csv")

# find the indices of the star1s and 2s that are in kepler-tgas
m = star1_df.source_id.isin(kplr_tgas.source_id).values * \
    star2_df.source_id.isin(kplr_tgas.source_id).values
star1_kic = star1_df.iloc[np.where(m)]  # list of star1s in kplr-tgas
star2_kic = star2_df.iloc[np.where(m)]  # list of star1s in kplr-tgas
print(np.shape(star1_kic))
print(star1_kic.keys())

plt.clf()
for i, _ in enumerate(star1_kic.ra):
    plt.plot([star1_kic.ra.values[i], star2_kic.ra.values[i]],
             [star1_kic.dec.values[i], star2_kic.dec.values[i]])
plt.plot(star1_kic.ra, star1_kic.dec, "r.")
plt.plot(star2_kic.ra, star2_kic.dec, "k.")
plt.savefig("kepler_ra_dec")

# do the same for the epic
m = star1_df.source_id.isin(epic_tgas.source_id).values * \
    star2_df.source_id.isin(epic_tgas.source_id).values
star1_epic = star1_df.iloc[np.where(m)]  # list of star1s in kplr-tgas
star2_epic = star2_df.iloc[np.where(m)]  # list of star1s in kplr-tgas
print(np.shape(star1_epic))

plt.clf()
plt.plot(star1_epic.ra, star1_epic.dec, "k.")
plt.plot(star2_epic.ra, star2_epic.dec, "k.")
plt.savefig("k2_ra_dec")
