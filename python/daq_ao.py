#!/usr/bin/env python2
import nidaqmx
task = nidaqmx.Task()
task.ao_channels.add_ao_voltage_chan('Dev1/ao0','mychannel',0,5)
task.start()
value = 1
task.write(value)
task.stop()
task.close()