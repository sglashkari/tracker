#!/usr/bin/env python2
import nidaqmx

SemiPeriod = nidaqmx.Task()
SemiPeriod.ci_semi_period_starting_edge = nidaqmx.constants.Edge.RISING
SemiPeriod.ci_channels.add_ci_semi_period_chan("Dev1/ctr0")
data_semiperiod = SemiPeriod.read(number_of_samples_per_channel=2)

duty_cycle = (data_semiperiod[0]/(data_semiperiod[0]+data_semiperiod[1]))*100
print("duty cycle = ", duty_cycle)