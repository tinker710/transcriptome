import os
import sys
import yaml

input_path=sys.argv[1]
output=sys.argv[2]
design_table=sys.argv[3]
gtf=sys.argv[4]
output_path=sys.argv[5]

with open(design_table,'r') as file1:
	data1=file1.readlines()[1:]

path1=output_path
file3=open(output,'w')
if not os.path.exists(path1):
	os.mkdir(path1)
for lines in data1:
	expt,ctrl=lines.strip().split("\t")
	path2=path1+expt+"_vs_"+ctrl+'/'
	os.mkdir(path2)
	with open(path2+'b1.txt','w') as file1:
		line1=input_path+expt+'_1.bam,'+input_path+expt+'_2.bam,'+input_path+expt+'_3.bam'
		file1.write(line1)
	with open(path2+'b2.txt','w') as file2:
		line2=input_path+ctrl+'_1.bam,'+input_path+ctrl+'_2.bam,'+input_path+ctrl+'_3.bam'
		file2.write(line2)
	line='nohup python /home/galaxy/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 '+path2+'b1.txt --b2 '+path2+'b2.txt --gtf '+gtf +' --od '+path2+" -t paired --readLength 100 --libType fr-firststrand &"
	file3.write(line)
	file3.write("\n")

file3.close()
print(1)


