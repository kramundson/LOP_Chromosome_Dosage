# USAGE: script parsed_mpileup_haploids.txt output.txt SNPfile.txt

import sys

mpup = open(sys.argv[1]) # mpup file
o = open(sys.argv[2]) # output file

## make a list of SNP position and their base call in the two parental genotypes

SNPfile = open(sys.argv[3]
SNPs = {}

SNPfile.readline() # headers

for mine in SNPfile:
    m = mine.split('\t')
    SNPChrom = m[0]
    SNPPos = m[1]
    SNPRef = m[2]
    SNP_A = m[3]
    SNP_B = m[4].spilt('\n')[0]
    SNPMega = int(SNPPos)/1000000
    if SNPChrom not in SNPs:
        SNPs[SNPChrom] = {}
    if SNPMega not in SNPs[SNPChrom]:
        SNPs[SNPChrom][SNPMega] = {}
    SNPs[SNPChrom][SNPMega][SNPPos] = [SNPRef,SNP_A,SNP_B]
SNPfile.close()

## Going through the mpileup now

indices = [] # need to keep indices, but not sure if need to keep the output file header
totalcount = 0 # cut
odd1 = 0 # cut 
odd2 = 0 # cut

line = mpup.readline() # header of mpup
l = line.split()
if line[0] == 'C':
    count = 0
    # header = 'Chrom\tPos\tRef\tA\tB' # cut
    header = '' # instead of getting chrom pos ref A B for file header, just print the indices
    for m in l:
        if m.split('-')[0] == 'Cov':
            indices.append(count)
            lib = m[4:]
            # header += '\tSnptype-'+lib+'\tSNP1-'+lib+'\tSNP2-'+lib+'\tTotalCov-'+lib+'\t%A-'+lib+'\tCovA-'+lib 
            header += lib+'\t' # header is now only the names of the libraries
        count += 1
    o.write(header+'\n')
    print indices
    
for line in mpup: # revise, keeping needed details
    l = line_split()
    chrom = l[0]
    pos = l[1]
    mega = int(pos)/1000000
    ref = l[2] # cut
    text = ''
    if chrom in SNPs:
        if mega in SNPs[chrom]:
            if pos in SNPs[chrom][mega]:
                a = SNPs[chrom][mega][pos][1]
                b = SNPs[chrom][mega][pos][2]
                refb = SNPs[chrom][mega][pos][0]
                if refb != ref:
                    print "ref problem"
                    break
                else:
                    # text = chrom+'\t'+str(pos)+'\t'+ref+'\t'+a+'\t'+b # cut
                    pass
                    
                for i in indices:
                    Cov = l[i]
                    if Cov == "."
                        text += '.\t'
                    else:
                        totalcount +=1
                        snp1 = l[i-3]
                        snp2 = l[i-2]
                        snp3 = l[i-1]
                        if snp3 != '.':
                            Snptype = 3
                        elif snp2 != '.':
                            Snptype = 2
                        else:
                            Snptype = 1
                        Snp1 = snp1.split('_')[0]
                        Cov1 = int(round(float(snp1.split("_")[1])*float(Cov)/100))
                        if snp2 == '.':
                            Snp2 = '.'
                        elif float(snp2.split('_')[1]) <= 10:
                            Snp2 = '.'
                        else:
                            Snp2 = snp2.split('_')[0]
                            Cov2 = int(round(float(snp2.split('_')[1])*float(Cov)/100))
                        if (Snp1 == a and Snp2 == b) or (Snp1 == b and Snp2 == a):
                            if Snp1 == a:
                                All = 0.5
                                CovA = Cov1
                            elif Snp1 == b:
                                All = 0.5
                                CovA = Cov2
                        elif Snp2 != '.':
                            All = '.'
                            odd2 += 1 # cut
                            CovA = '.'
                        elif Snp1 == b:
                            All = 0
                            CovA - 0
                        elif Snp1 == a:
                            All = 1
                            CovA = Cov
                        else:
                            All = '.'
                            odd1 += 1
                            CovA = '.'
                        text += chrom+'_'+pos+'\t' # THIS IS THE MOST IMPORTANT PART
                text += '\n'
                o.write(text)
o.close()
mpup.close()

print totalcount
print odd1
print odd2