#!/usr/bin/env python2
import nidaqmx
import time
task_write = nidaqmx.Task()
task_write.ao_channels.add_ao_voltage_chan('Dev1/ao0','mychannel',0,5)
task_write.start()
#task_read = nidaqmx.Task()
#task_read.ai_channels.add_ai_voltage_chan("Dev1/ai0")
#task_read.start()
start=0; stop=6; increment=1
for k in range(start, stop, increment):
	value = k
	if value>5:
		value=5
	task_write.write(value)
	time.sleep(1)
#	value = task_read.read()
#	print(round(value,2))
task_write.stop()
task_write.close()
#task_read.stop()
#task_read.close()