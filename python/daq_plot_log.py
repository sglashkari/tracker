#!/usr/bin/env python2
import nidaqmx
import time
import numpy as np
import matplotlib.pyplot as plt

# Initialize Logging
Tstop = 2 # Logging Time [seconds]
Ts = 0.001 # Sampling Time [seconds]
N = int(Tstop/Ts)
data = []

# Initialize DAQ Device
task = nidaqmx.Task()
task.ai_channels.add_ai_voltage_chan("Dev1/ai0")
task.start()

# Open File
file = open("tempdata.txt", "w")

# Write Data Function
def writefiledata(t, x):
	time = str(t)
	value = str(round(x, 2))
	file.write(time + "\t" + value)
	file.write("\n")

# Logging Temperature Data from DAQ Device
for k in range(N):
	value = task.read()
	print(round(value,1))
	data.append(value)
	time.sleep(Ts)
	writefiledata(k*Ts, value)
	print(k*Ts)

# Terminate DAQ Device
task.stop()
task.close()

# Close File
file.close()

# Plotting
t = np.arange(0,Tstop,Ts)
plt.plot(t,data, "-o")
plt.xlabel('t [s]')
plt.grid()
plt.show()
#block=False)
#plt.pause(3)
#plt.close()
