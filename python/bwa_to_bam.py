import sys
import optparse
from multiprocessing import Pool
from subprocess import Popen, PIPE, call

def runcommand(cmd):
    print >>sys.stderr, "calling "+cmd
    try:
        retcode = call(cmd, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
    except OSError as e:
        sys.exit ("Execution of "+cmd+" failed: "+str(e))


def runscript(sample_string):
    if sample_string.strip() == "":
        return
    else:
        print >>sys.stdout, "sample", sample_string
        host,sample,location = sample_string.split()

    if lines != 0:
        p1 = Popen(["samtools", "view", location], stdout=PIPE, stderr=logfile)
        p2 = Popen(["head", "-n", str(lines)], stdin=p1.stdout, stdout=PIPE)

        smallbamfilename = str(sample+".small.bam")
        smallbamfile = open(smallbamfilename, "w")
        p3 = Popen(["samtools", "view", "-S", "-u", "-"], stdin=p2.stdout, stdout=smallbamfile, stderr=logfile)
        p3.communicate()
        smallbamfile.close()
        p2.terminate()
        runbwa(smallbamfile, sample)
        return
    else:
        runbwa(location,sample)

def runbwa(bamfile,sample):
        cmd = "bwa aln -b1 %s %s > %s.1.sai" % (refname,smallbamfilename,sample)
        print >>sys.stdout, cmd
        runcommand(cmd)
        cmd = "bwa aln -b2 %s %s > %s.2.sai" % (refname,smallbamfilename,sample)
        runcommand(cmd)
        cmd = "bwa sampe %s %s.1.sai %s.2.sai %s.small.bam %s.small.bam > %s.sam" % (refname,sample,sample,sample,sample,sample)
        runcommand(cmd)
        cmd = "rm %s.1.sai; rm %s.2.sai; rm %s.small.bam" % (sample,sample,sample)
        runcommand(cmd)
        cmd = "samtools view -S -b -u -o %s.bam %s.sam" % (sample,sample)
        runcommand(cmd)
        # remove unmapped pairs
        cmd = "samtools view -F 4 -b %s.bam > %s.reduced.bam" % (sample,sample)
        runcommand(cmd)
        cmd = "rm %s.sam" % (sample)
        runcommand(cmd)
        cmd = "mv %s.reduced.bam %s.bam" % (sample,sample)
        runcommand(cmd)
        # sort pairs
        cmd = "samtools sort %s.bam %s.sorted" % (sample,sample)
        runcommand(cmd)
        cmd = "rm %s.bam" % (sample)
        runcommand(cmd)


global refname
global lines
#Parse Command Line
parser = optparse.OptionParser()
parser.add_option("-i", "--input", type="string", default="", dest="input", help="A list of files to run script on")
parser.add_option("-r", "--reference", type="string", default="~/Populus/reference_seqs/populus.trichocarpa.cp.fasta", dest="ref", help="The reference genome")
parser.add_option("-p", "--processes", default=2, type="int", dest="processes", help="Number of processes to use")
parser.add_option("-n", "--number", default=5000, type="int", dest="num", help="Number of short reads to use")

(options, args) = parser.parse_args()
refname = options.ref
lines = options.num

if options.input == "":
    sys.exit("Sample file must be provided.\n")

global logfile
logfile = open (str(options.input+".log"), "w")

try:
    open(options.ref, "r").close()
    cmd = "bwa index %s" % (refname)
    runcommand(cmd)
except IOError as e:
    sys.exit("Reference file " + options.ref + " not found\n")

print >>sys.stdout, "using "+refname+" with "+str(lines)

pool = Pool(int(options.processes))

#read the location file
try:
    handle = open(options.input, "r")
    samples = []
    for line in handle:
        sample = line.rstrip()
        samples.append(sample)
    handle.close()
except IOError as e:
    sys.exit("Sample file " + options.input + " not found\n")

pool.map(runscript, samples)

logfile.close()
