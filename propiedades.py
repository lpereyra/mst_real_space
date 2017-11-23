import numpy as np
import matplotlib.pyplot as plt
import sys

bintype = np.dtype([
        ("flag", np.int32),
        ("size", np.int32),
        ("len", np.float32),
        ("elong", np.float32),
        ("rms", np.float32),
                  ])

binario = open(sys.argv[1],"rb")
N = np.fromfile(binario, dtype=np.int32, count=1)
N = N[0]
print "Num Segmentos",N

data = np.fromfile(binario, dtype=bintype)
data["rms"] /= data["len"]

num_bins = 100
mask = (data["flag"]>=0)
print "Mask Segmentos",sum(mask)
mask = (data["size"]>3)*mask
print "Segmentos mask N>2",len(data["size"][mask])

#HISTOGRAMA LONGITUD
log = 0
len_min = (data["len"][mask].min())-1e-10
len_max = (data["len"][mask].max())+1e-10
print "len min", len_min
print "len max", len_max
if log==0:
  len_bins = np.linspace(len_min,len_max,num=num_bins)
else:
  len_bins = np.logspace(np.log10(len_min),np.log10(len_max),num=num_bins)

#HISTOGRAMA ELONGACION
log = 0;
elong_min = (data["elong"][mask].min())-1e-10
elong_max = (data["elong"][mask].max())+1e-10
print "elong min", elong_min
print "elong max", elong_max
if log==0:
  elong_bins = np.linspace(elong_min,elong_max,num=num_bins)
else:
  elong_bins = np.logspace(np.log10(elong_min),np.log10(elong_max),num=num_bins)

#HISTOGRAMA RMS
log = 0;
rms_min = (data["rms"][mask].min())-1e-10
rms_max = (data["rms"][mask].max())+1e-10
print "rms min", rms_min
print "rms max", rms_max
if log==0:
  rms_bins = np.linspace(rms_min,rms_max,num=num_bins)
else:
  rms_bins = np.logspace(np.log10(rms_min),np.log10(rms_max),num=num_bins)

f, axarr = plt.subplots(3)

#HISTOGRAMA LONGITUD
axarr[0].hist(data["len"][mask],len_bins,color="r",normed=True)
axarr[0].set_xlabel("Longitud")
axarr[0].set_xlim([len_min,len_max])

#HISTOGRAMA ELONGACION
axarr[1].hist(data["elong"][mask],elong_bins,color="g")
axarr[1].set_xlabel("Elongacion")
axarr[1].set_xlim([elong_min,elong_max])

#HISTOGRAMA RMS
axarr[2].hist(data["rms"][mask],rms_bins,color="b")
axarr[2].set_xlabel("RMS")
axarr[2].set_xlim([rms_min,rms_max])

f.subplots_adjust(hspace=0.5)

plt.show()


