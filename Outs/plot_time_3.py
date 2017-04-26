import matplotlib.pyplot as plt
import numpy as np

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = float(rect.get_height())
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%.2f' % float(height),
                ha='center', va='bottom')


cores = [68]
arr= [20]
rng=[100.0]
dims=[]
for z in range(16,257,16):
	dims.append(z)
elems=[65536]
resSeq=[]
res68=[]
res68e=[]
with open('serialnovec_3.out') as file:#with open('serialnovec.out') as file:
	list1 = file.readlines()
with open('testnovec_3.out') as file:
	list2 = file.readlines()
with open('final_vec_mem_3.out') as file:
	list3 = file.readlines()
for i in range(0,len(dims)):
	f0=0.0
	tmp = int((list1[i].split()[0]).split('=')[1])
	#for ctr in range(0,5):
	f0+= float((list1[i].split()[3]).split('=')[1])
	f0=f0/tmp
	resSeq.append(f0)
	f1=0.0
	tmp = int((list2[5*i].split()[0]).split('=')[1])
	for ctr in range(0,5):
		f1+= float((list2[5*i+ctr].split()[3]).split('=')[1])
	res68.append((f1/5.0/tmp))
	f1=0.0
	tmp = int((list3[5*i].split()[0]).split('=')[1])
	for ctr in range(0,5):
		f1+= float((list3[5*i+ctr].split()[3]).split('=')[1])
	res68e.append((f1/5.0/tmp))
fig, ax = plt.subplots(figsize=(16, 9),dpi=100)
N = 16
ind = np.arange(N)  # the x locations for the groups
width = 0.4      # the width of the bars	
#rectsSeq = ax.bar(ind , resSeq, width, color='r',label='Seq')
#rectsSeq = ax.bar(ind , resSeq, width, color='r',label='Sequential')
rects68 = ax.bar(ind, res68, width, color='r',label='68 cores baseline')
rects68e = ax.bar(ind + width +0.05 , res68e, width, color='y',label='68 cores final')
ax.set_ylabel('Time/Iter (sec/iter)')
ax.set_xlabel('Dimensions')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(dims)
#autolabel(rects68)
#autolabel(rects68e)
#ax.set_ylim([0,max(res1+res2+res3+res4+res5+res6+res7+res8+res9+res10+res11)*5/3])
ax.set_title("Varied dimension number")
ax.legend(prop={'size':10},loc=2)
plt.savefig("Varied_dimension_number.png", bbox_inches='tight')
plt.close()

