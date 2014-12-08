#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@package    Find_overlap_read
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# STANDARD IMPORTS
import os
from time import time
import optparse
import csv

# THIRD PARTY IMPORTS
#try:
    #import pysam # from pysam 0.8.1

#except ImportError as E:
    #print (E)
    #print ("Please verify your dependencies. See Readme for more informations\n")
    #exit()

# PACKAGE IMPORTS
from Interval import Interval

#~~~~~~~MAIN FUNCTION~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
    """

    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self):

        # Define instance variables
        self.program_name =  "find_overlap_read"
        self.program_version = "0.1"

        # Parse CLI arguments
        self.opt_dict = self._optparser()
        for key, value in self.opt_dict.items():
            print (key, "=", value)

        # Parse the csv file containing interval coordinates
        with open(self.opt_dict["interval_file"], newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:

                # Verifying if the row fields are valid and skip the row if errors were found
                # If valid = create a class tracked instance of Interval
                try:
                    if len(row) < 3:
                        raise ValueError ("Not enough values in the row")
                    elif len(row) == 3:
                        Interval(ref_name=row[0], start=int(row[1]), end=int(row[2]))
                    else:
                        Interval(ref_name=row[0], start=int(row[1]), end=int(row[2]), name=row[3])

                except ValueError as E:
                    print (E, "\tSkiping row")

        Interval.printInstances()
        self.interval_list = Interval.getInstances()

    def __str__(self):
        msg = "MAIN CLASS\n"
        msg+= "\tParameters list\n"
        for i, j in self.__dict__.items():
            msg+="\t{}\t{}\n".format(i, j)
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        """
        stime = time()
        return (stime)

    ##~~~~~~~PRIVATE METHODS~~~~~~~#

    def _optparser(self):
        """
        Parse command line arguments with optparse and verify the file validity
        @param program_name Name of the program
        @param program_version Version of the program
        @return A dictionnary containing:
        """
        # Usage and version
        usage_string = ("{} -f genomic_interval.csv [-o Output_prefix] [-b/-r] f1.bam(sam),"
        "[f2.bam(sam)...fn.bam(sam)\n Parse a BAM/SAM file(s) and extract reads overlapping given"
        " genomic coordinates\n").format(self.program_name)
        version_string = "{} {}".format(self.program_name, self.program_version)
        optparser = optparse.OptionParser(usage = usage_string, version = version_string)

        # Define optparser options
        optparser.add_option('-f', dest="interval_file",
        help="Path of the tab separated file contaning genomic interval (mandatory)")
        optparser.add_option( '-o', '--output', default="out", dest="output_prefix",
        help="Facultative option to indicate the name of the output prefix (default = out)" )
        optparser.add_option('-b', '--no_bam', action="store_false", dest="bam_output", default=True,
        help="Don't output bam file(s) (default = True)")
        optparser.add_option('-r', '--no_report', action="store_false", dest="report_output", default=True,
        help="Don't output report file(s) (default = True)")

        # Parse arg and return a dictionnary_like object of options
        options, args = optparser.parse_args()

        try:
            # Verify if mandatory opt were given
            if not args:
                raise ValueError ("No alignment file (bam/sam) was provided")
            if not options.interval_file:
                raise ValueError ("No alignment file (bam/sam) was provided")

            # Verify files readability
            for fp in [options.interval_file]+args :
                if not os.access(fp, os.R_OK):
                    raise ValueError ("{} is not a valid file".format(fp))

        except ValueError as E:
            print (E)
            optparser.print_help()
            exit (1)

        # Create a dictionnary of options for further ease to use
        return ({'interval_file' : options.interval_file,
                'output_prefix' : options.output_prefix,
                'bam_output' : options.bam_output,
                'report_output' : options.report_output,
                'align_files' : args })


    #def _sam_spliter (self):
        #"""
        #"""
        #with pysam.Samfile(self.sam, "r" ) as samfile:
            #self.bam_header = samfile.header

            ## Give the header of the sam file to all Reference.Instances to respect the same order
            ## references in sorted bam files
            #Reference.set_global("bam_header", self.bam_header)

            ## Create a dict to collect unmapped and low quality reads
            #Secondary = Sequence (name = 'Secondary', length = 0)
            #Unmapped = Sequence (name = 'Unmapped', length = 0)
            #LowMapq = Sequence (name = 'LowMapq', length = 0)
            #self.garbage_read = [Secondary, Unmapped, LowMapq]

            #for read in samfile:
                ## Always remove secondary alignments
                #if read.is_secondary:
                    #Secondary.add_read(read)
                ## Filter Unmapped reads
                #elif read.tid == -1:
                    #Unmapped.add_read(read)
                ## Filter Low MAPQ reads
                #elif read.mapq < self.min_mapq:
                    #LowMapq.add_read(read)
                ## Finally if all is fine attribute the read to a Reference
                #else:
                    #Reference.addRead(samfile.getrname(read.tid), read)

        ## Removing the original sam file which is no longer needed
        #remove(self.sam)
        #self.sam = None






#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    main = Main()
    main()
