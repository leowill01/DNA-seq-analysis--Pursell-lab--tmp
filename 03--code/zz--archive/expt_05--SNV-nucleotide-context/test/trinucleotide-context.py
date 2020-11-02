# NOTE: Usage: $ python trinucMutContext.py [VCF_FILE] [REF_GENOME_CHR_FASTAS_DIRECTORY] > newVCF.vcf.trint

# Written by Nate Ungerleider

from Bio import SeqIO
import sys
import os
import glob

def rcomp(seq):

    dna = {'A': 'T', 'C': 'G', 'T': 'A', 'G':'C'}
    return ''.join(reversed([dna[i.upper()] if i.upper() in dna else 'X' for i in seq]))

'''

Takes a vcf file (chromosome, position, id, ref, alt, etc in that order) and directory of fastas as input. Prints vcf file with added column for trinucleotide context to standard out.

'''

def parse():

    script_name = os.path.basename(sys.argv.pop(0))

    helpmsg = "\nUsage: {script} vcf_file genome_directory\n".format(script=script_name)
    helpmsg += "\ngenome_directory must contain individual fasta files for each chromosome."
    helpmsg += " i.e. 'chr1.fa', 'chr2.fa', etc. Or '1.fa', '2.fa'; as long as the text up until the '.fa' matches"
    helpmsg += " the chromosome names in the vcf file.\n"


    if len(sys.argv) < 2:
        print("\nYou need to supply at least two arguments..a vcf file and a directory (containing genome fasta files.)\n")
        sys.exit(helpmsg)

    if sys.argv[1] in ['-h', '--help']:
        sys.exit(helpmsg)

    vcf = sys.argv.pop(0)
    if not os.path.exists(vcf):
        sys.exit("VCF file not found: %s" % vcf)

    genome_dir = sys.argv.pop(0)

    fastas = glob.glob(os.path.join(genome_dir,'*fa'))
    if not fastas:
        sys.exit("No fasta files found in your genome directory: %s" % genome_dir)

    return vcf, genome_dir


def process_vcf(vcf, genome_dir):

    c_init = ''
    with open(vcf, 'r') as fh:
        for line in fh:
            line = line.strip('\n')

            if line[0] == '#':
                print(line)
                continue

            chromosome, position, _, ref, alt, *_ = line.split('\t', maxsplit=5)
            position = int(position)
            if chromosome != c_init:
                c_init = chromosome
                try:
                    loaded_fa = SeqIO.read(os.path.join(genome_dir, c_init + '.fa'), 'fasta')
                except:
                    sys.exit("%s not found or not in the proper format" % c_init)
            position = int(position)
            if len(ref) != 1 or len(alt) != 1:
                print(line)
                continue
            subseq = str(loaded_fa[position - 2: position + 1].seq).upper()
            if ref in 'GA':
                ref = rcomp(ref)
                alt = rcomp(alt)
                subseq = rcomp(subseq)
            mutation = "{ref}>{alt}-{subseq}".format(ref=ref, alt=alt, subseq=subseq)

            line += '\t%s' % mutation
            print(line)


vcf, genome_dir = parse()
process_vcf(vcf, genome_dir)
