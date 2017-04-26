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
dims=[16]
rng=[100.0]
arr=[]
for z in range(10,211,20):
	arr.append(z)
elems=[65536]
resSeq=[]
res68=[]
res68e=[]
with open('serialnovec_2.out') as file:#with open('serialnovec.out') as file:
	list1 = file.readlines()
with open('testnovec_2.out') as file:
	list2 = file.readlines()
with open('final_vec_mem_2.out') as file:
	list3 = file.readlines()
for i in range(0,len(arr)):
	f0=0.0
	tmp = int((list1[i].split()[0]).split('=')[1])
	#for ctr in range(0,5):
	f0+= float((list1[i].split()[3]).split('=')[1])
	f0 =f0/tmp
	resSeq.append(f0)
	f1=0.0
	tmp = int((list2[5*i].split()[0]).split('=')[1])
	for ctr in range(0,5):
		f1+= float((list2[5*i+ctr].split()[3]).split('=')[1])
	res68.append(f0/(f1/5.0/tmp))
	f1=0.0
	tmp = int((list3[5*i].split()[0]).split('=')[1])
	for ctr in range(0,5):
		f1+= float((list3[5*i+ctr].split()[3]).split('=')[1])
	res68e.append(f0/(f1/5.0/tmp))
fig, ax = plt.subplots(figsize=(16, 9),dpi=100)
N = 11
ind = np.arange(N)  # the x locations for the groups
width = 0.4      # the width of the bars	
#rectsSeq = ax.bar(ind , resSeq, width, color='r',label='Seq')
#rectsSeq = ax.bar(ind , resSeq, width, color='r',label='Sequential')
rects68 = ax.bar(ind, res68, width, color='g',label='68 cores baseline')
rects68e = ax.bar(ind + width +0.05 , res68e, width, color='b',label='68 cores final')
ax.set_ylabel('Speedup')
ax.set_xlabel('Cluster Centers')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(arr)
#autolabel(rects68)
#autolabel(rects68e)
#ax.set_ylim([0,max(res1+res2+res3+res4+res5+res6+res7+res8+res9+res10+res11)*5/3])
ax.set_title("Varied Cluster number Speedup(4MB)")
ax.legend(prop={'size':10},loc=2)
plt.savefig("Varied_cluster_number_speedup(4MB).png", bbox_inches='tight')
plt.close()

