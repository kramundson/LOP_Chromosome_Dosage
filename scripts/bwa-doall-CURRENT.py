#! /usr/bin/env python
#! /usr/bin/env python
import sys, math, os, time
from optparse import OptionParser
from collections import defaultdict
import re, subprocess

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2017
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#This program processes library FASTQ files through .sorted.bam files.
#This process includes the following:
#1. If option is selected to remove chimeric reads, this is done first. New library files are generated
#   with nc appended for no chimeric, the orignial fastq files are still kept, and a statistic file
#   is generated for the removal called rescan-cut-log.txt
#2a. If a paired end overamp is specified, the INTERLEAVED pair ended Fastq file is then split to forward/reverse
#2b. Specified overamp (single or paired) is performed after bwa alignment is done SINGLE ENDED.
#3. Overamp is run, this gives unique and aligned numbers in the file master-OverAmp.txt.
#   If the -o/--overamp option is used, the unique reads sam files will be used to continue instead of the original
#4. The .bam files are generated, .sorted.bam are made from these. The orignial unsorted bams are removed
#
#INPUT:
#This program must be run in a folder of .fq fastq files. It will look at the contents of the folder and run on all
#files that end with ".fq".
#
#OUTPUT:
#Folders for sam (if used), non-overamp sams, sai, sorted bams, bai, original fq files, uncut fq files (if remove chimeric is used) 
#are generated and relevent files moved to the correct result directory.
#
#
#The sorted.bam files that are generated here can then be used to create a mpileup file (using our mplieup script package), 
#which is then used for mutation detection / genotyping (using our MAPS package).





from optparse import OptionParser
usage = "USAGE: bwa-samtools.py -d database_file.fa [-c] [-p] [-o] [-t #threads] <Note overamp can not be down with bwa meme or bwa sw>"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--database-file", dest="database", help="Input database file for mapping.")
parser.add_option("-c", "--chimeric", dest="chimeric",  action="store_true", default = False, help="Remove Chimeric Reads")
parser.add_option("-o", "--overamp", dest="overamp",  action="store_true", default = False, help="Remove Overamplified Reads")
#parser.add_option("-O", "--supressoveramp", dest="nooveramp",  action="store_false", default = True, help="Suppress all overamp actions, default = True")
parser.add_option("-m", "--mode", dest="mode", type = "str", default='s', help="When mapping as SE (no insert estimation), Read type, s = SE (default), p = PE no salvage, ps = PE with SE salvage (keep the one ampped read), pb = PE with PE salvage(keep both reads only if one maps)")
parser.add_option("-D", "--disableSD", dest="disableSD", action="store_false", default = True, help="If running in ps or pb mode, this turns off collecting pairs mapping in the same direction (0,0) or (16,16)")
parser.add_option("-t", "--thread", dest="threads", default="8", help="How many threads to use during alignment.")
parser.add_option("-q", "--trimqual", dest="trimqual", default="20", help="Default mapping quality, for use as the bwa aln -q X during alignment.")
parser.add_option("--bwa", "-b", dest="pathBWA",  type = "str", default='', help="File path to BWA")
parser.add_option("--samtools", "-a", dest="pathSAM",  type = "str", default='', help="File path to Samtools")
parser.add_option("--scripts", "-s", dest="pathScript",  type = "str", default='/home/mclieberman/projects/web-pipelines/bwa-package/', help="File path to all package scripts. (if using overamp!)")
parser.add_option("-B", "--backtrack", dest="memmode", action="store_false", default = True, help="Default: Use the bwa (backtrack) algorithm instead of mem defualt (bwa-mem) ")
parser.add_option("-S", "--sw", dest="swmode", action="store_true", default = False, help="Use the bwa sw algorithm instead of mem default")
parser.add_option("-P", "--pair", dest="pairmode", action="store_true", default = False, help="Use true pair ended mapping, (insert size calculated in mapping) must be sued with bwa mem")
parser.add_option("-X", "--XoutModule", dest="modules", action="store_false", default = True, help="Loading modules for a unix module absed system, defualt = False ")


(opt, args) = parser.parse_args()



if opt.modules == False:
   if not os.environ.has_key('MODULE_VERSION'):
   	os.environ['MODULE_VERSION_STACK'] = '3.2.10'
   	os.environ['MODULE_VERSION'] = '3.2.10'
   else:
   	os.environ['MODULE_VERSION_STACK'] = os.environ['MODULE_VERSION']
   os.environ['MODULESHOME'] = '/software/modules/3.2.10/x86_64-linux-ubuntu14.04/Modules/3.2.10'
   
   if not os.environ.has_key('MODULEPATH'):
   	f = open(os.environ['MODULESHOME'] + "/init/.modulespath", "r")
   	path = []
   	for line in f.readlines():
   		line = re.sub("#.*$", '', line)
   		if line is not '':
   			path.append(line)
   	os.environ['MODULEPATH'] = ':'.join(path)
   
   if not os.environ.has_key('LOADEDMODULES'):
   	os.environ['LOADEDMODULES'] = ''
   	
   def module(*args):
   	if type(args[0]) == type([]):
   		args = args[0]
   	else:
   		args = list(args)
   	(output, error) = subprocess.Popen(['/software/modules/3.2.10/x86_64-linux-ubuntu14.04/Modules/%s/bin/modulecmd' % os.environ['MODULE_VERSION'], 'python'] + 
   			args, stdout=subprocess.PIPE).communicate()
   	exec output
   	
   module("load samtools/0.1.19")
   #module("load bwa/0.7.13")



#if opt.swmode == True:
#   nooveramp = True


opt.mode = str.lower(opt.mode)
if opt.mode not in ['s', 'p', 'ps', 'pb'] and opt.pairmode == False:
   parser.error("Please check your command line mode paramter, must be s, p, ps, or pb")
   

if opt.disableSD == True:
   DISaddon = ''
else:
   DISaddon = ' -D '


database = opt.database
remchim = opt.chimeric
threads = opt.threads
trimQual = opt.trimqual
pathBWA = opt.pathBWA
pathSam = opt.pathSAM
pathElse = opt.pathScript

#from Isabelle Henry, this looks for chimeric reads to be removed
def rescanCutter(fname):
   print "Removing chimeric"
   all = os.listdir(os.getcwd())
   fqs = filter(lambda x: ".fq" in x or ".fastq" in x, all)
   fqs.sort()
   out = open(fname,'w')
   out.write('File\tIntact\tCut\tRemoved\tTotal\t%intact\t%cut\t%removed\n')
   os.system("mkdir orig-fq")
   for fq in fqs:
      print str(fq)
      f = open(fq)
      newname = str(fq).split('.')[0]+'nc.fq'
      o = open(str(fq).split('.')[0]+'nc.fq','w')
      countyes = 0
      countrem = 0
      countno = 0 
      while True:
         name = f.readline()
         if name == "":
            break
         seq = f.readline()
         plus = f.readline()
         qual = f.readline()
         if 'CATG' in seq[4:]:
            seq = seq[:4]+seq[4:].split('CATG')[0]+'CATG'+'\n'
            if len(seq) > 35:
               o.write(name+seq+plus+qual[:len(seq)-1]+'\n')
               countyes +=1
            else:
               countrem +=1
         else:
            o.write(name+seq+plus+qual)
            countno +=1
      f.close()
      o.close()
   
      total = countno + countyes + countrem
      data = [fq, countno, countyes, countrem, total, round(100.00*countno/total,3), round(100.00*countyes/total,3), round(100.00*countrem/total,3)]
      data = map(lambda x: str(x), data)
      out.write('\t'.join(data)+'\n')   
      os.system("mv "+fq+" orig-fq/")
   out.close()

def uninterleave(fname):
   namet = ''.join(fname.split('.')[:-1])
   u = open(fname)
   f = open(namet+"-1.fq",'w')
   r = open(namet+"-2.fq",'w')
   while True:
      name1 = u.readline()
      if name1 == '':
         break
      seq1 = u.readline()
      p1 = u.readline()
      qual1 = u.readline()
      name2 = u.readline()
      seq2 = u.readline()
      p2 = u.readline()
      qual2 = u.readline()
      f.write(name1+seq1+p1+qual1)
      r.write(name2+seq2+p2+qual2)   
   f.close()
   r.close()
   u.close()

#/////////// do index for bwa if have not
if '/' in database:
   dbfiles = os.listdir(database[:database.rindex('/')+1])
   dname = database[database.rindex('/')+1:]
else:
   dbfiles = os.listdir(os.getcwd())
   dname = database
ind = filter(lambda x: dname in x, dbfiles)

index = ""
dsize = os.path.getsize(database)
if dsize < 1048576000:
   index = " is "
else:
   index = " bwtsw "   
#makes sure ref has been indexed, will skip if already done
ref = ['', '.pac', '.ann', '.amb', '.bwt', '.sa']
check = list((x in ref for x in map(lambda x: x[len(database):],ind)))
if check.count('False') > 1 or len(check) < len(ref):
   print(pathBWA+"bwa index -a"+index+" "+database) 
   os.system(pathBWA+"bwa index -a"+index+" "+database)


#/////////////////////////

if remchim == True:
   rescanCutter("rescan-cut-log.txt")

#setup result directories
li = os.listdir(os.getcwd())
if opt.memmode == True or opt.swmode == True or opt.pairmode == True:
   dirs = ['sam', 'bam', 'bai', 'fq', 'usam']
   if False in list((x in li for x in dirs)):
      os.system("mkdir sam bam bai fq usam")
else:
   dirs = ['sam', 'sai', 'bam', 'bai', 'fq', 'usam']
   if False in list((x in li for x in dirs)):
      os.system("mkdir sam sai bam bai fq usam")



if opt.overamp == False:
   os.system("rm -r usam")



#make header for overamp run
if opt.overamp == True:
   temp = open('master-OverAmp.txt', 'w')
   temp.write('\t'.join(["File", "Type", "#Non-Clonal", "#Aligned", "%NonClonal", "All", "%Mapped"])+'\n')
   temp.close()
else:
   temp = open("master-counts.txt", 'w')
   temp.write('\t'.join(["File", "Type", "#Aligned", "All", "%Mapped"])+'\n')
   temp.close()   

todo = filter(lambda x: x[-3:] == '.fq' or x.endswith("fastq"), li)
todo.sort()
print "Aligning"
for file in todo:
   print file
   name = file.split('.')[0]
   if opt.pairmode == True:
      os.system(pathBWA+"bwa mem -t "+threads+" -p "+database+" "+file+" > "+name+"_aln.sam")
      print(pathBWA+"bwa mem -t "+threads+" -p "+database+" "+file+" > "+name+"_aln.sam")
      if opt.overamp == True:
         print(pathElse+"overamp-6.py -f "+name+"_aln.sam -o u"+name+"_aln.sam "+ " -m tp >> master-OverAmp.txt")
         os.system(pathElse+"overamp-6.py -f "+name+"_aln.sam -o u"+name+"_aln.sam "+ " -m tp >> master-OverAmp.txt")
         name = 'u'+name
      else:
         f = open(name+"_aln.sam")
         o = open("master-counts.txt",'a')
         flagcount = defaultdict(int)
         pereads = 0
         pealign = 0
         while 1:
            x = f.readline()
            if x == "":
               break
            if x[0] == '@':
               continue
            l = x.split('\t')
            #other alignment
            if int(l[1]) > 200:
               #print x
               continue               
            y = f.readline()
            m = y.split('\t')
            while int(m[1]) > 200:
               y = f.readline()
               m = y.split('\t')      
            
            flagset =[l[1],m[1]]
            if flagset == ['99','147'] or flagset == ['83','163']:
               pereads+=1
               pealign+=1
            elif flagset == ['147','99']:
               x[333333]
            else:
               #print flagset
               pereads+=1     
      
         o.write('\t'.join([name, "TP", str(pealign), str(pereads), str(float(pealign)/pereads*100)])+'\n')
         o.close()   
         
   else:
      if opt.memmode == False and opt.swmode == False:
         os.system(pathBWA+"bwa aln -t "+threads+" -q "+trimQual+" "+database+" "+file+" > "+name+"_aln.sai")
         os.system(pathBWA+"bwa samse "+database+" "+name+"_aln.sai "+file+" > "+name+"_aln.sam")
      elif opt.memmode == True:
         os.system(pathBWA+"bwa mem -t "+threads+" "+database+" "+file+" > "+name+"_aln.sam")
      else:
         os.system(pathBWA+"bwa bwasw -t "+threads+" "+database+" "+file+" > "+name+"_aln.sam") 
      if opt.overamp == True:
         print(pathElse+"overamp-6.py -f "+name+"_aln.sam -o u"+name+"_aln.sam "+ " -m "+opt.mode+DISaddon+" >> master-OverAmp.txt")
         os.system(pathElse+"overamp-6.py -f "+name+"_aln.sam -o u"+name+"_aln.sam "+ " -m "+opt.mode+DISaddon+" >> master-OverAmp.txt")
      else:
         f = open(name+"_aln.sam")
         o = open("master-counts.txt",'a')
         flagcount = defaultdict(int)
         sereads = 0
         pereads = 0
         sealign = 0
         pealign = 0
         while 1:
            x = f.readline()
            if x == "":
               break
            if x[0] == '@':
               continue
            l = x.split('\t')
            if opt.mode == 's':
               if l[1] in ['0', '16']:
                  sealign+=1
                  sereads+=1
               else:
                  sereads+=1                       
            else:
               #other alignment
               if int(l[1]) > 200:
                  #print x
                  continue               
               y = f.readline()
               m = y.split('\t')
               while int(m[1]) > 200:
                  y = f.readline()
                  m = y.split('\t')
                  
               flagset =[l[1],m[1]]
               if flagset == ['0','16'] or flagset == ['16','0']:
                  pereads+=1
                  pealign+=1
               elif flagset == ['0','0'] or flagset == ['16','16']:
                  if opt.disableSD == False and opt.mode in ['ps','pb']:
                     pereads+=1
                     pealign+=1
                  else:
                     pereads+=1
               elif '0' in flagset or '16' in flagset:
                  if opt.mode in ['ps','pb']:
                     sereads+=1
                     sealign+=1
               else:
                  pereads+=1
                  if flagset != ['4','4']:
                     x[43456789765456]
               
         if opt.mode == 's':   
            o.write('\t'.join([name, opt.mode, str(sealign), str(sereads), str(float(sealign)/sereads*100)])+'\n')
         elif opt.mode == 'p':
            o.write('\t'.join([name, opt.mode, str(pealign), str(pereads), str(float(pealign)/pereads*100)])+'\n')
         elif opt.mode in ['ps','pb']:
            o.write('\t'.join([name, opt.mode+'-SE', str(sealign), str(sereads), str(float(sealign)/sereads*100)])+'\n')    
            o.write('\t'.join([name, opt.mode+'-PE', str(pealign), str(pereads), str(float(pealign)/pereads*100)])+'\n')
         o.close()
         
         
         
         
         
         
         
         
         
      if opt.overamp == True:
         name = 'u'+name


   os.system(pathSam+"samtools view -bS "+name+"_aln.sam"+" > "+name+"_aln.bam")
   os.system(pathSam+"samtools sort "+name+"_aln.bam"+" "+name+"_aln.sorted")
   os.system("rm -f "+name+"_aln.bam")
   os.system(pathSam+"samtools index "+name+"_aln.sorted.bam")
   if opt.overamp == True:
      os.system('mv u*_aln.sam usam/') 
   os.system('mv *_aln.sam sam/')     
   #os.system("mv "+file+" fq/") 
   if opt.memmode == False and opt.swmode == False and opt.pairmode == False:
      os.system("mv *.sai sai/") 
   os.system("mv *.bam bam/")
   os.system("mv *.bai bai/")








