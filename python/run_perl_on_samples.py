import sys, os
import optparse
from socket import gethostname
from multiprocessing import Pool

scriptname = "";

def assign_script( script ):
	global scriptname
	scriptname = script;
	print("scriptname is %s\n" % scriptname);

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def runscript(sample_string):
    cmd = "perl %s -input %s" % (scriptname,sample_string)
    print(cmd)
    os.system(cmd)

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", default=None, dest="input",
                      help="A list of files to run script on")
    parser.add_option("-s", "--script", default=None, dest="script",
                      help="The script to batch")
    parser.add_option("-p", "--processes", default=None, dest="processes",
                      help="Number of processes to use")
    (options, args) = parser.parse_args()

    try:
        open(options.input, "r").close()
    except TypeError, e:
        stop_err("You need to supply the input file:\n" + str(e))

    try:
        open(options.script, "r").close()
    except TypeError, e:
        stop_err("Script file not found:\n" + str(e))

    print ("options.script\n")
    assign_script(options.script)

    pool = Pool(processes=int(options.processes))

    #read the location file
    handle = open(options.input, "r")
    samples = []
    for line in handle:
        sample = line.rstrip()
        samples.append(sample)
    handle.close()

    #this can be run in parallel
    pool.map(runscript, samples)

if __name__=="__main__":
    __main__()
