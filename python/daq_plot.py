#!/usr/bin/env python2
import nidaqmx
import time
import numpy as np
import matplotlib.pyplot as plt
# Initialize Logging
Tstop = 10 # Logging Time [seconds]
Ts = 1 # Sampling Time [seconds]
N = int(Tstop/Ts)
data = []
# Initialize DAQ Device
task = nidaqmx.Task()
task.ai_channels.add_ai_voltage_chan("Dev1/ai0")
task.start()
#Logging Temperature Data from DAQ Device
for k in range(N):
	value = task.read()
	print(round(value,1))
	data.append(value)
	time.sleep(Ts)
# Terminate DAQ Device
task.stop()
task.close()
# Plotting
t = np.arange(0,Tstop,Ts)
plt.plot(t,data, "-o")
plt.title('Voltage')
plt.xlabel('t [s]')
plt.ylabel('Voltage [V]')
plt.grid()
#Tmin = 18; Tmax = 28
plt.show()