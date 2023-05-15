#!/usr/bin/env python2
import nidaqmx
task = nidaqmx.Task()
task.ai_channels.add_ai_voltage_chan("Dev1/ai0")
task.start()
value = task.read()
print(value)
task.stop
task.close()