import os
import sys
import yaml

input_path=sys.argv[1]
output=sys.argv[2]
config_file=open(sys.argv[3])
config=yaml.load(config_gile)
a1=config["design_table1"]

ps.mkdir(input_path)
os.chdir(input_path)
rmats=open(output,'w')
with open(a1,'r') as file1:
	data1=file1.readlines()[1:]
for lines in data1:
	expt, ctrl=lines.strip().split("\t")
	path2=expt+"_vs_"+ctrl+"/"
	os.mkdir(input_path+path2+)
	with open(input_path+path2+'b1.txt') as file1:
		line='bam/'+expt+'_1.bam,bam/'+expt+'_2.bam,bam/'+expt+'_3.bam'
		file1.write(line)
	with open(input_path+path2+'b2.txt') as file2:
		line='bam/'+ctrl+'_1.bam,bam/'+ctrl+'_2.bam,bam/'+ctrl+'_3.bam'
		file2.write(line)
	line='python /home/galaxy/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 '+input_path+path2+'\b1.txt --b2 '+input_path+path2+'\b2.txt --gtf '+config['gtf']+'-od '+input_path+path2+" -t paired +readlength 100 --libTypr fr-firststrand"


