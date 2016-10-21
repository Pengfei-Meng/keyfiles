
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

import pickle
import pdb
import numpy as np


gxs = [64, 128, 256, 512, 1024]
nlevel = 4

nl_iter = []
prc_calls = []
cost_solve = []
nl_iter_np = []
prc_calls_np = []
cost_solve_np = []



fig = plt.figure()
clr = ['bv-', 'mo-', 'gx-', 'cd-', 'r*-']
clr_np = ['bv-.', 'mo-.', 'gx-.', 'cd-.', 'r*-.']
labels = ['64','128','256','512', '1024']
i = 0


#///////////////////////////////////
with open('64WtP', 'rb') as f:
	whole_64 = pickle.load(f)
error_norm_64 = whole_64['error_norm']	
del error_norm_64[0]
error_norm_64 = error_norm_64/error_norm_64[0]
f.close()

with open('128WtP', 'rb') as f:
	whole_128 = pickle.load(f)
error_norm_128 = whole_128['error_norm']	
del error_norm_128[0]
error_norm_128 = error_norm_128/error_norm_128[0]
f.close()

with open('256WtP', 'rb') as f:
	whole_256 = pickle.load(f)
error_norm_256 = whole_256['error_norm']	
del error_norm_256[0]
error_norm_256 = error_norm_256/error_norm_256[0]
f.close()

with open('512WtP', 'rb') as f:
	whole_512 = pickle.load(f)
error_norm_512 = whole_512['error_norm']	
del error_norm_512[0]
error_norm_512 = error_norm_512/error_norm_512[0]
f.close()

with open('1024WtP', 'rb') as f:
	whole_1024 = pickle.load(f)
error_norm_1024 = whole_1024['error_norm']	
del error_norm_1024[0]
error_norm_1024 = error_norm_1024/error_norm_1024[0]
f.close()

#//////////////////////////////////

with open('64NoP', 'rb') as f:
	whole_64np = pickle.load(f)
error_norm_64np = whole_64np['error_norm']	
del error_norm_64np[0]
error_norm_64np = error_norm_64np/error_norm_64np[0]
f.close()


with open('128NoP', 'rb') as f:
	whole_128np = pickle.load(f)
error_norm_128np = whole_128np['error_norm']	
del error_norm_128np[0]
error_norm_128np = error_norm_128np/error_norm_128np[0]
f.close()

with open('256NoP', 'rb') as f:
	whole_256np = pickle.load(f)
error_norm_256np = whole_256np['error_norm']	
del error_norm_256np[0]
error_norm_256np = error_norm_256np/error_norm_256np[0]
f.close()

with open('512NoP', 'rb') as f:
	whole_512np = pickle.load(f)
error_norm_512np = whole_512np['error_norm']	
del error_norm_512np[0]
error_norm_512np = error_norm_512np/error_norm_512np[0]
pdb.set_trace()
f.close()

with open('1024NoP', 'rb') as f:
	whole_1024np = pickle.load(f)
error_norm_1024np = whole_1024np['error_norm']	
del error_norm_1024np[0]
error_norm_1024np = error_norm_1024np/error_norm_1024np[0]
f.close()



fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
pl64, = plt.plot(error_norm_64, clr[0])
pl128, = plt.plot(error_norm_128, clr[1])
pl256, = plt.plot(error_norm_256, clr[2])
pl512, = plt.plot(error_norm_512, clr[3])
pl1024, = plt.plot(error_norm_1024, clr[4])


pl64_np, = plt.plot(error_norm_64np, clr_np[0])
pl128_np, = plt.plot(error_norm_128np, clr_np[1])
pl256_np, = plt.plot(error_norm_256np, clr_np[2])
pl512_np, = plt.plot(error_norm_512np, clr_np[3])
pl1024_np, = plt.plot(error_norm_1024np, clr_np[4])


plt.xlabel('Nonlinear Iteration', fontsize = 14)
plt.ylabel('Normalized Objective Function Value', fontsize = 14)
plt.legend([pl64, pl128, pl256, pl512, pl1024, pl64_np, pl128_np, pl256_np, pl512_np, pl1024_np],   #  pl512_np  pl1024_np   '(512,512, No P)', (1024,1024, No P)',
	['(64,64)','(128,128)','(256,256)','(512,512)', '(1024,1024)','(64,64, No P)','(128,128, No P)','(256,256, No P)', '(512,512, No P)', '(1024,1024, No P)'], loc =1, prop={'size':8})
# plt.legend([pl64, pl128, pl256, pl512, pl1024], ['64','128','256','512','1024'], loc =3, prop={'size':14})
ax1.set_yscale('log')
plt.savefig('sum_converg_temp.eps')

i += 1

#-----------------------

for gx in gxs:

	gy = gx 
	caseid = str(gx) + "WtP"       # + "_" + str(gy) + "_" + str(nlevel) 

	with open(caseid, 'rb') as f: 
		whole = pickle.load(f)

	error_norm = whole['error_norm']	
	del error_norm[0]

	nl_iter.append(len(error_norm)-1)

	prec_n = whole['precond_calls'] - whole['setup_precond']
	prc_calls.append(prec_n)

	#-----------------------

	time_system = whole['time_solve_system']
	del time_system[0]

	cost = (whole['kona_time'] - whole['setup_time'])/ np.mean(time_system)

	cost_solve.append(cost)

	f.close()

##############-------------------------
gxsp = [64, 128, 256, 512, 1024]
for gx in gxsp:

	gy = gx 
	caseid = str(gx)+"NoP"    # + "_" + str(gy) + "_" + str(nlevel) + "noP"

	with open(caseid, 'rb') as f:
		whole = pickle.load(f)

	error_norm = whole['error_norm']	
	del error_norm[0]

	nl_iter_np.append(len(error_norm)-1)

	prec_n = whole['precond_calls'] - whole['setup_precond']
	prc_calls_np.append(prec_n)

	#-----------------------

	time_system = whole['time_solve_system']
	del time_system[0]

	cost = (whole['kona_time'] - whole['setup_time'])/ np.mean(time_system)

	cost_solve_np.append(cost)

	f.close()


#######################################

clr2 = ['b', 'm', 'g', 'c', 'r']
mk = [u'v', u'o', u'x', u'd', u'*']
# s = [20*2**n for n in range(len(gxs))];

fig2 = plt.figure()
fig2.subplots_adjust(bottom=0.2)
ax2 = fig2.add_subplot(1,1,1)
# plt.plot(gxs, nl_iter, 'bo-', gxs, nl_iter_np,  'bo-.', linewidth = 1.5)
ni_, = plt.plot(gxs, nl_iter, 'bo-', linewidth = 1.5)
ni_np, = plt.plot(gxsp, nl_iter_np,  'bs-', linewidth = 1.5)
plt.legend([ni_, ni_np],['Preconded RSNK', 'not Preconded RSNK'],loc =1, prop={'size':14})
# plt.ylim(ymin = 0, ymax=10)
plt.xticks(gxs, labels, rotation=0, fontsize = 14)
plt.ylabel('Number of NonLinear Iterations', fontsize = 14)
plt.xlabel('Number of Design Variables', fontsize = 14)
# plt.legend()
plt.savefig('sum_NLIter_temp.eps')


fig3 = plt.figure()
fig3.subplots_adjust(bottom=0.2)
ax3 = fig3.add_subplot(1,1,1)
# plt.plot(gxs, prc_calls, 'mv-', gxs, prc_calls_np, 'mv-.', linewidth = 1.5)
prc_, = plt.plot(gxs, prc_calls, 'bo-', linewidth = 1.5)
prc_np, = plt.plot(gxsp, prc_calls_np, 'bs-', linewidth = 1.5)
plt.legend([prc_, prc_np],['Preconded RSNK', 'not Preconded RSNK'], loc =1, prop={'size':14})
# plt.ylim(ymin = 0, ymax=550)
plt.xticks(gxs, labels, rotation=0, fontsize = 14)
plt.ylabel('No of Preconditioner calls', fontsize = 14)
plt.xlabel('Number of Design Variables', fontsize = 14)
plt.savefig('sum_PrecondNo_temp.eps')

fig4 = plt.figure()
fig4.subplots_adjust(bottom=0.2)
ax4 = fig4.add_subplot(1,1,1)
# plt.plot(gxs, cost_solve, 'g*-', gxs, cost_solve_np, 'g*-.', linewidth = 1.5)
ct_, = plt.plot(gxs, cost_solve, 'bo-', linewidth = 1.5)
ct_np, = plt.plot(gxsp, cost_solve_np, 'bs-', linewidth = 1.5)
plt.legend([ct_, ct_np], ['Preconded RSNK', 'Not Preconded RSNK'], loc =1, prop={'size':14} )
# plt.ylim(ymin = 0, ymax = 100)
plt.xticks(gxs, labels, rotation=0, fontsize = 14)
plt.ylabel('Cost - Number of function analysis', fontsize = 14)
plt.xlabel('Number of Design Variables', fontsize = 14)
plt.savefig('sum_cost_temp.eps')









