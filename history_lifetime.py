import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse
from itertools import groupby
from datetime import datetime

import lifetimeMeasurement
from selection import track_dtype

parser = argparse.ArgumentParser()
parser.add_argument('infileList',
                    nargs = 1,
                    type = str,
                    help = 'flattened hit data from selection.py')
parser.add_argument('--plotfile', '-p',
                    default = "",
                    type = str,
                    help = 'optional file to which resulting lifetime history plot is saved')
args = parser.parse_args()

hitData, trackData = lifetimeMeasurement.read_files(args.infileList)
# sort for groupby
trackData = np.sort(trackData, order=["unixtime"])

# Divide tracks into time bins
bin_hours = 4 #bin width of history [hours], should be a devisor of 24

lifetimes = []
datetimes = []
for _, this_group in  groupby(trackData, 
                    key = lambda x: (datetime.fromtimestamp(x["unixtime"]).month,
                                     datetime.fromtimestamp(x["unixtime"]).day,
                                     datetime.fromtimestamp(x["unixtime"]).hour // bin_hours)):
    trackInfo = np.array(list(this_group), dtype=track_dtype)
    lifetimeMeasurement.calc_dQdx(hitData, trackInfo)

    if len(lifetimeMeasurement.driftTime) > 100:
        # fit the observation to a decay model
        fitModel = lifetimeMeasurement.expModelWithNorm((100, 100))
        fitModel.fit((lifetimeMeasurement.driftTime, lifetimeMeasurement.dQdxs))
        time_start = min(trackInfo["unixtime"])
        time_end = max(trackInfo["unixtime"])
        datetimes.append(datetime.fromtimestamp((time_start+time_end)/2))
        lifetimes.append(fitModel.bf_tau())

plt.rcParams["font.size"] = 15
plt.plot(datetimes, lifetimes, marker='x')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d\n%H:%M'))
plt.xlabel('Date')
plt.ylabel(r'Lifetime [$\mu$s]')
plt.ylim(bottom=0)
plt.tight_layout()
plt.grid(True)
if args.plotfile:
    plt.savefig(args.plotfile)
else:
    plt.show()