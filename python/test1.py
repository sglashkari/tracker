#!/usr/bin/env python2
import nidaqmx

device = nidaqmx.system.device.Device("Dev1")
for tr in device.terminals:
	print(tr)