import numpy as np
import struct

with open("data.data",mode="rb") as file:
    data=file.read()

struct.unpack("iiiii",data[:20])
data=struct.unpack("i" * ((len(data) -24) // 4), data[20:-4])
print data

"""
dt=np.dtype([('t',float)])
data=np.fromfile("data.dat",dtype=dt)

print data.shape
"""
