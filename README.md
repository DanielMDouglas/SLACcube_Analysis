# SLACcube_Analysis

Simple selection and analysis scripts for analyzing the data collected by the SLAC single cube prototype system.

## Selection

Selection is currently very simple and consists of a veto region or "margin" around the outside of the anode plane.  If the beginning of a given track intersects with this region, it is likely that it is in fact a wall-piercing track and its absolute drift position is uncertain.

The selection script can be run on a larpix file processed throught the h5flow pipeline.  It may be called like:

```
python selection.py [-h] [-p PLOTFILE] -o OUTFILE infileList [infileList ...]
```

This will compile all of the given input files into a single output file.  It may be advantageous to do this with only one input file at a time.

The output format is a simplified flat array stored in an HDF5 file.  The 3D position, estimated charge, and an ad-hoc track index is stored in "hits" object, and the track informations like total charge, colinear parameter (= width / length), length, cosine of the polar angle measured from drift axis, unixtime, and the track index corresponding to the one in "hits" are stored in "track" object.

The script also plots basic features of tracks. It will be saved in a file if a file name is specified by `-p` option

## Lifetime Measurement

The electron lifetime can be measured using one of these files, as demonstrated in the `lifetimeMeasurement.py` script.  It may be called like:

```
python lifetimeMeasurement.py [-h] [--plotfile PLOTFILE] infileList
```

This script analyses the relationship between the position along the drift axis and the dQ/dx (Q is represented by ADC counts) of small segments along tracks to examine the extinction rate of charge as it drifts.  This data is compared to an exponential decay model with a normalization term.  The result is displayed directly using matplotlib, or saved to a file whose location is provided by the user. Basic features (dQ, dx, colinear parameter (= width / length), and dQ/dx) for each segment are also shown if the resultant plot is displayed directly. Calibration constant from ADC counts to incident particle energy [MeV] is also calculated.

Also, lifetime values for binned periods can be plotted with `history_lifetime.py` like:

```
python history_lifetime.py [-h] [--plotfile PLOTFILE] infileList
```

History of electron lifetime value is plotted directly in a window, or saved in a file if specified by user.

## Track Display

Also included as an example is a simple track display that can be used to parse selected events by eye.  Its usage is:

```
python trackDisplay.py [-h] [-t TRACKID] [-v] [--plotfile PLOTFILE] infile
```

The user provides a file and a trackID to display and the track is shown along with an outline the the detector volume directly using matplotlib, or saved to a file whose location is provided by the user. Track information is provided when `-v` option is specified.

To see many tracks continuously, type:

```
python continuousDisplay.py [-h] [-t TRACKID] [-v] [--plotfile PLOTFILE] infile
```

A track is displayed like trackDisplay.py. Then type 'k' or 'n' for the next track, 'j' or 'b' for the previous track, or 'e' or 'q' to quit.

## Contributing

Contributions are welcome!  Please feel free to fork and submit pull requests, issues, or reach out to me directly via email or Slack.  Chances are, if you've found this page, you know how to find me!  Happy analyzing!
