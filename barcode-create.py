# name = input('Enter file name: ')
import math
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plot_url
import sys

bzeroes = []
dzeroes = []
bzeroesinf = []
dzeroesinf = []
bones = []
dones = []
bonesinf = []
donesinf = []
btwoers = []
dtwoers = []
btwoersinf = []
dtwoersinf = []

x = []
maxx = [-1,-1,-1]

if len(sys.argv) == 1:
    print 'need filename'
    sys.exit(-1)
# file = input('Enter filename','s')
print(sys.argv[1])
with open(sys.argv[1], 'r') as file:
	for line in file:
		l = line.split(" ")
		if(l[2][:-1]=='inf'):
			l[2] = -(10)
		l[0] = float(l[0])
		l[1] = float(l[1])
		l[2] = float(l[2])
		# print((l[0]),(l[1]),l[2])

		if(l[0]==0):	# zero dim homology
			if(l[2]<0):
				bzeroesinf.append(l[1])
				dzeroesinf.append(l[2])
			else:#(bzeroes == [] or (bzeroes[-1] != (l[2]-l[1]))):	#check if same value is entered
				bzeroes.append(l[1])
				dzeroes.append(l[2])
				if(l[2]-l[1]>maxx[0]):
					maxx[0] = l[2]-l[1]
		elif (l[0]==1):
			if(l[2]<0):
				bonesinf.append(l[1])
				donesinf.append(l[2])
			else:#if(bones == [] or (bones[-1] != (l[2]-l[1]))):
				bones.append(l[1])
				dones.append(l[2])
				if(l[2]-l[1]>maxx[1]):
					maxx[1] = l[2]-l[1]
		else:
			if(l[2]<0):
				btwoersinf.append(l[1])
				dtwoersinf.append(l[2])
			else:#(btwoers == [] or (btwoers[-1] != (l[2]-l[1]))):
				btwoers.append(l[1])
				dtwoers.append(l[2])
				if(l[2]-l[1]>maxx[2]):
					maxx[2] = l[2]-l[1]

# print(twoers)
for ind, val in enumerate(dzeroesinf):
	dzeroesinf[ind] = maxx[0]+10
for ind, val in enumerate(donesinf):
	donesinf[ind] = maxx[1]+10
for ind, val in enumerate(dtwoersinf):
	dtwoersinf[ind] = maxx[2]+10


# delete noise based on input parameter
todelete = len(bzeroes)+len(bzeroesinf)- int(sys.argv[2])
print todelete, len(bzeroesinf)
if todelete>0:
	del bzeroes[0:todelete]
	del dzeroes[0:todelete]
# delete noise based on input parameter
todelete = len(bones)+len(bonesinf)- int(sys.argv[2])
print todelete
if todelete>0:
	del bones[0:todelete]
	del dones[0:todelete]
# delete noise based on input parameter
todelete = len(btwoers)+len(btwoersinf) - int(sys.argv[2])
print todelete
if todelete>0:
	del btwoers[0:todelete]
	del dtwoers[0:todelete]


z = [len(bzeroes),len(bones),len(btwoers)]
bzeroes.extend(bzeroesinf)
bones.extend(bonesinf)
btwoers.extend(btwoersinf)
dzeroes.extend(dzeroesinf)
dones.extend(donesinf)
dtwoers.extend(dtwoersinf)

bzeroesinf = bzeroes
bonesinf = bones
btwoersinf = btwoers
dzeroesinf = dzeroes
donesinf = dones
dtwoersinf = dtwoers

dzeroesinf = list(reversed(dzeroesinf))
bzeroesinf = list(reversed(bzeroesinf))
y_pos = np.arange(len(dzeroesinf))
barlist0 = []
barlist1 = []
barlist2 = []
plt.figure(1)
plt.subplot(211)
plt.xlabel('0-dim')
for element in y_pos:
	if element>= len(bzeroesinf) - z[0]:
		barlist0.append(plt.plot([bzeroesinf[element], dzeroesinf[element] ], [element , element], color='red', linestyle='-', linewidth=4))
	else:
		barlist0.append(plt.plot([bzeroesinf[element], dzeroesinf[element] ], [element , element], color='black', linestyle='-', linewidth=4))

bonesinf = list(reversed(bonesinf))
donesinf = list(reversed(donesinf))
y_pos = np.arange(len(bonesinf))
plt.subplot(212)
plt.xlabel('1-dim')
for element in y_pos:
	if element>= len(bonesinf) - z[1]:
		barlist1.append(plt.plot([bonesinf[element], donesinf[element]], [element, element], color='blue', linestyle='-', linewidth=4))
	else:
		barlist1.append(plt.plot([bonesinf[element], donesinf[element]], [element, element], color='black', linestyle='-', linewidth=4))

plt.show()
