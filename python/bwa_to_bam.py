import sys, os
import optparse
from socket import gethostname
from multiprocessing import Pool

refname = "~/Populus/reference_seqs/populus.trichocarpa.cp.fasta";

def assign_ref( ref ):
	global refname
	refname = ref;
	print("refname is %s\n" % refname);

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def runscript(sample_string):
    """
    remote filters the bam file
    """
    host,sample,location = sample_string.split()

    #This happens locally
    cmd = "bwa aln %s -b1 %s > %s.1.sai" % (refname,location,sample)
    sys.stderr.write( "%s\n" % cmd )
    os.system(cmd)
    cmd = "bwa aln %s -b2 %s > %s.2.sai" % (refname,location,sample)
    os.system(cmd)
    cmd = "bwa sampe %s %s.1.sai %s.2.sai %s %s > %s.sam" % (refname,sample,sample,location,location,sample)
    os.system(cmd)
    cmd = "samtools view -S -b -u -o %s.bam %s.sam" % (sample,sample)
    os.system(cmd)
    cmd = "rm %s.1.sai; rm %s.2.sai; rm %s.sam" % (sample,sample,sample)
    os.system(cmd)

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", default=None, dest="input",
                      help="A list of files to run script on")
    parser.add_option("-r", "--reference", default=None, dest="ref",
                      help="The reference genome")
    parser.add_option("-p", "--processes", default=None, dest="processes",
                      help="Number of processes to use")
    (options, args) = parser.parse_args()

    try:
        open(options.input, "r").close()
    except TypeError, e:
        stop_err("You need to supply the input file:\n" + str(e))

    try:
        open(options.ref, "r").close()
    except TypeError, e:
        stop_err("Reference file not found:\n" + str(e))

    print ("options.ref\n")
    assign_ref(options.ref)

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
