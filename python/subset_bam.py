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

        p1 = Popen(["samtools", "view", location], stdout=PIPE, stderr=logfile)
        p2 = Popen(["head", "-n", str(lines)], stdin=p1.stdout, stdout=PIPE)

        smallbamfilename = str(sample+".small.bam")
        smallbamfile = open(smallbamfilename, "w")
        p3 = Popen(["samtools", "view", "-S", "-u", "-"], stdin=p2.stdout, stdout=smallbamfile, stderr=logfile)
        p3.communicate()
        smallbamfile.close()
        p2.terminate()

global lines
#Parse Command Line
parser = optparse.OptionParser()
parser.add_option("-i", "--input", type="string", default="", dest="input", help="A list of files to run script on")
parser.add_option("-p", "--processes", default=2, type="int", dest="processes", help="Number of processes to use")
parser.add_option("-n", "--number", default=5000, type="int", dest="num", help="GB to subset")

(options, args) = parser.parse_args()
lines = options.num * 500000

if options.input == "":
    sys.exit("Sample file must be provided.\n")

global logfile
logfile = open (str(options.input+".log"), "w")

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
