import numpy as np


reported = list()
with open('report_quantiles.txt') as ins:
    for line in ins:
        if len(line)>3 and line[0]!="%":
            ls = line.split()
            reported.append(map(float,ls[2:]))



with open('single_attached.txt') as ins:
    data = list()
    frame = list()
    for line in ins:
        if "% end" in line:
            data.append(frame)
            frame = list()
        elif len(line) > 3 and line[0] != "%":

            ls = line.split()
            frame.append(float(ls[7]))

quantile_num = [0.25,.5,.75]

for frame in [-1,-10,-20]:
    print "---- frame %u ----" % frame
    for i,q in enumerate(quantile_num):

        quantile_value = np.quantile(data[frame],q)
        abscissas = data[frame]
        abscissas = np.sort(abscissas)
        min_id= np.argmin(np.abs(abscissas-quantile_value))
        print abscissas[min_id], reported[frame][i]

