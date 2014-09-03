import sys, os
import optparse
from socket import gethostname
from multiprocessing import Pool

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def runscript(sample_string):
    host,sample,location = sample_string.split()
    pathname = os.path.dirname(sys.argv[0])
    head,tail = os.path.split(pathname)
    cmd = "bash %s/converting/bam_to_fasta.sh %s %s" % (head, location, sample)
    os.system(cmd)

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", default=None, dest="input",
                      help="A list of files to run script on")
    parser.add_option("-p", "--processes", default=1, dest="processes",
                      help="Number of processes to use")
    (options, args) = parser.parse_args()

    try:
        open(options.input, "r").close()
    except TypeError, e:
        stop_err("You need to supply the input file:\n" + str(e))

    pool = Pool(processes=int(options.processes))

    #read the location file
    handle = open(options.input, "r")
    samples = []
    for line in handle:
        sample = line.rstrip()
        samples.append(sample)
    handle.close()

    pool.map(runscript, samples)

if __name__=="__main__":
    __main__()

