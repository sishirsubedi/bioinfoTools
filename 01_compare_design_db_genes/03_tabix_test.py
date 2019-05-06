from subprocess import Popen, PIPE


def bgzip(filename):
    """Call bgzip to compress a file."""
    Popen(['bgzip', '-f', filename])

def tabix_index(filename,
        preset="gff", chrom=1, start=4, end=5, skip=0, comment="#"):
    """Call tabix to create an index for a bgzip-compressed file."""
    Popen(['tabix', '-p', preset, '-s', chrom, '-b', start, '-e', end,
        '-S', skip, '-c', comment])

def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    for line in process.stdout:
        print( line.strip().split())

# tabix_query("clinvar_mod.vcf.gz",'chr1',1,955555)
tabix_query("/home/hhadmin/gnomad.exomes.r2.0.2.sites.vcf.gz",'chr1',1,1055555)
