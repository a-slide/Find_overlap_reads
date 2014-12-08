#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@package    Find_overlap_read
@brief      Main file of the program, containing a Main 
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
try:
    import pysam # from pysam 0.8.1

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()

# PACKAGE IMPORTS
from Interval import Interval


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
    """
    
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, interval_file=None, align_files=[], output_prefix="out",
        bam_output=True, report_output=True):

        # Define instance variables
        self.program_name =  "find_overlap_read"
        self.program_version = "0.1"
        self.max_tlen = 1000

        # For interactive IDE call = simple parsing without verification of files 
        if interval_file:
            self.opt_dict = {
                'interval_file' : interval_file,
                'output_prefix' : output_prefix,
                'bam_output' : bam_output,
                'report_output' : report_output,
                'align_files' : align_files }
            
        # Else for Shell call = parse and verify more throrougly CLI arguments
        else:
            self.opt_dict = self._optparser()
        
        # Print the dictionnary
        for key, value in self.opt_dict.items():
            print (key, "=", value)

        # Parse the csv file containing interval coordinates
        print ("Parsing the CSV csv file containing interval coordinates")
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
                
                # Can by catched if not enough values or invalid numeric value in row
                except ValueError as E:
                    print (E, "\tSkiping row")
        
        if Interval.countInstances() == 0:
            raise ValueError ("No valid row found in the interval file. Please see readme file")
        
        Interval.printInstances()
        print ("{} valid interval(s) found in {}".format(
            Interval.countInstances(),
            self.opt_dict["interval_file"]))

        self.interval_list = Interval.getInstances()
        #Interval.resetInstances()

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
        
        """"""
        stime = time()
        
        # Iterate over the bam/sam alignment files given as positional arguments 
        for align_file in self.opt_dict["align_files"]:
            sam = pysam.AlignmentFile(align_file, "rb")
            
            tot_read = map_read = 0
            
            # Iterate over the aligned reads in sam
            for read in sam:
                tot_read+=1
                
                # If the read is mapped only
                if not read.is_unmapped:
                    read_ref = sam.getrname(read.reference_id)
                    map_read+=1
                    
                    # Iterate over all instances of interval
                    for interval in self.interval_list:
                        
                        # If the read overlap the coordinates, add to the interval
                        if interval.is_overlapping (read_ref, read.reference_start,
                        read.reference_end):
                            interval.add_read(read)
                        
                        # Else try to see if the pair overlap the coordinates
                        elif not read.mate_is_unmapped:
                            mate_ref = sam.getrname(read.next_reference_id)
                            
                            # Verify first if the pair of read is a properly mapped pair
                            if self.pair_is_concordant(
                                rref=read_ref,
                                mref=mate_ref,
                                tlen=read.template_length,
                                ror=read.is_reverse,
                                mor=read.mate_is_reverse):
                                
                                # If the pair overlap the coordinates, add to the interval
                                if interval.is_overlapping(read_ref, read.reference_start,
                                read.next_reference_start):
                                    interval.add_read(read)
                        
            if self.opt_dict["report_output"]:
                self.write_report(tot_read, map_read, Interval.get_report())
                
            if self.opt_dict["bam_output"]:
                self.write_bam(sam.header, Interval.get_read())
            
            Interval.printInstances()
            Interval.resetReadCount()
            Interval.resetReadList()
            
            #print ("Total reads :", tot_read)
            #print ("Mapped reads :", map_read)
            #print ("Overlapping reads :", overlap_read)
            #print ("Overlapping pair :", overlap_pair)
            
        return (stime)
        
        
        
    def pair_is_concordant(self, rref, mref, tlen, ror, mor):
        """
        """
        # Verify lenght of the template, if read and mate are on the same reference and if
        # their orientation is concordant (R1F2 or F1R2) 
        
        return 0 < abs(tlen) <= self.max_tlen and rref == mref and ror != mor

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


#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    main = Main()
    main()
