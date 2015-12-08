The directory mf2005 is meant to put MODFLOW examples that were explicit for MF2005.

However, most problemsn can be solved with eigther mf2k or mf2005. However mf2005 lacks
the parameter optimization and calibration facilities of mf2k. On the other hand mf2005
does not have an artificial limit of the maximum number of stress periods in a simulation.
This is important for the exmaple Dutchtop, where 9 years is simulated on a daily basis (3288)
stress periods to extract the mean yearly maximum and minimum groundwater elevations, which are
computed over 8 years of data.

Note that you an easily switch between mf2k and mf205 by setting the MF5 flag in the NAM worksheet
accompanying the Dutchtop model.


TO 100823