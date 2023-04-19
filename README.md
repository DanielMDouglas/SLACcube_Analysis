# SLACcube_Analysis

Simple selection and analysis scripts for analyzing the data collected by the SLAC single cube prototype system.

## Selection

Selection is currently very simple and consists of a veto region or "margin" around the outside of the anode plane.  If the beginning of a given track intersects with this region, it is likely that it is in fact a wall-piercing track and its absolute drift position is uncertain.

The selection script can be run on a larpix file processed throught the h5flow pipeline.  It may be called like:

```
python selection.py [-h] -o OUTFILE infileList [infileList ...]
```

This will compile all of the given input files into a single output file.  It may be advantageous to do this with only one input file at a time.

The output format is a simplified flat array stored in an HDF5 file.  Only the 3D position, estimated charge, and an ad-hoc track index is stored.

## Lifetime Measurement

The electron lifetime can be measured using one of these files, as demonstrated in the `lifetimeMeasurement.py` script.  It may be called like:

```
./lifetimeMeasurement.py [-h] [--plotfile PLOTFILE] [--use_absolute] infileList [infileList ...]
```

This script analyses the relationship between the position along the drift axis and the charge (relative to the initial hit or in absolute terms) to examine the extinction rate of charge as it drifts.  This data is compared to an exponential decay model (with or without a normalization term, depending upon the charge metric).  The result is displayed directly using matplotlib, or saved to a file whose location is provided by the user.

## Track Display

Also included as an example is a simple track display that can be used to parse selected events by eye.  Its usage is:

```
./trackDisplay [-h] [-t TRACKID] [--plotfile PLOTFILE] infile
```

The user provides a file and a trackID to display and the track is shown along with an outline the the detector volume directly using matplotlib, or saved to a file whose location is provided by the user.

## Contributing

Contributions are welcome!  Please feel free to fork and submit pull requests, issues, or reach out to me directly via email or Slack.  Chances are, if you've found this page, you know how to find me!  Happy analyzing!
