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
from Find_overlap_read_src.Interval import Interval
from pyDNA3.Utilities import file_basename, file_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
    """
    
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, interval_file=None, template_len=1000, align_files=[],
        bam_output=True, report_output=True):

        # Define instance variables
        self.program_name =  "find_overlap_read"
        self.program_version = "0.1"
        self.max_tlen = 1000

        # For interactive IDE call = simple parsing without verification of files 
        if interval_file:
            self.opt_dict = {
                'interval_file' : interval_file,
                'template_len' : template_len,
                'bam_output' : bam_output,
                'report_output' : report_output,
                'align_files' : align_files }
            
        # Else for Shell call = parse and verify more throrougly CLI arguments
        else:
            self.opt_dict = self._optparser()

        # Parse the csv file containing interval coordinates
        print ("Parsing the file containing genomic interval coordinates\n")
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
            raise ValueError ("No valid row found in the interval file. Please see readme file\n")
        
        #Interval.printInstances()
        print ("{} valid interval(s) found in {}\n".format(
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
        """
        """
        stime = time()
        
        # Iterate over the bam/sam alignment files given as positional arguments 
        for al_file in self.opt_dict["align_files"]:
            print ("Analyse file {}".format(file_name(al_file)))
            sam = pysam.AlignmentFile(al_file, "rb")
            tot_read = map_read = 0
            report_list = []
            
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
                            if self._pair_is_concordant(
                                read_ref,
                                mate_ref,
                                read.is_reverse,
                                read.mate_is_reverse,
                                read.reference_start,
                                read.next_reference_start):
                                
                                # If the pair overlap the coordinates, add to the interval
                                if interval.is_overlapping(read_ref, read.reference_start,
                                read.next_reference_start):
                                    interval.add_read(read)
            
            if self.opt_dict["report_output"]:
                print("Write report")
                self._write_report(al_file, tot_read, map_read, Interval.get_report())
                
            if self.opt_dict["bam_output"]:
                ("Write Bam")
                self._write_bam(al_file, sam.header, Interval.get_read())
            
            Interval.printInstances()
            Interval.resetReadCount()
            Interval.resetReadList()
            
        print ("\n##### DONE #####\n")
        print ("Total execution time = {}s".format(round(time()-stime, 2)))
        

    ##~~~~~~~PRIVATE METHODS~~~~~~~#

    def _optparser(self):
        """
        Parse command line arguments with optparse and verify the file validity
        @param program_name Name of the program
        @param program_version Version of the program
        @return A dictionnary containing:
        """
        # Usage and version
        usage_string = ("{} -f genomic_interval.csv [-b/-r] f1.bam(sam),"
        "[f2.bam(sam)...fn.bam(sam)\n Parse a BAM/SAM file(s) and extract reads overlapping given"
        " genomic coordinates\n").format(self.program_name)
        version_string = "{} {}".format(self.program_name, self.program_version)
        optparser = optparse.OptionParser(usage = usage_string, version = version_string)

        # Define optparser options
        optparser.add_option('-f', dest="interval_file",
        help="Path of the tab separated file contaning genomic interval (mandatory)")
        optparser.add_option('-t', dest="template_len", default=1000,
        help="Maximum length of the template between 2 paired reads to be considered as "
        "a concordant pair (default = 1000)")
        optparser.add_option('-b', '--no_bam', action="store_false", dest="bam_output",
        default=True, help="Don't output bam file(s) (default = True)")
        optparser.add_option('-r', '--no_report', action="store_false", dest="report_output",
        default=True, help="Don't output report file(s) (default = True)")
        

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
                 'template_len' : int(options.template_len),
                 'bam_output' : options.bam_output,
                 'report_output' : options.report_output,
                 'align_files' : args })

    def _pair_is_concordant(self, read_ref, mate_ref, read_reverse, mate_reverse, read_start, mate_start):
        """
        Verify lenght of the template, if read and mate are on the same reference and if
        their orientation is concordant (R1F2 or F1R2)
        """
        if read_ref == mate_ref:
            if not read_reverse and mate_reverse:
                if 0 < mate_start-read_start <= self.opt_dict['template_len']:
                    return True
            elif read_reverse and not mate_reverse:
                if 0 < read_start-mate_start <= self.opt_dict['template_len']:
                    return True
                    
        # In all other conditions
        return False
        
    def _write_report(self, al_file, tot_read, map_read, interval_report):
        """
        Write a brief report containing the number of read mapped for each interval
        """
        # Generate a output name for the report
        outname = "./{}_report.txt".format(file_basename(al_file))
        with open(outname, "w") as report:
            report.write("{}\n".format(file_basename(al_file)))
            report.write("Total_reads\t{}\n".format(tot_read))
            report.write("Mapped_reads\t{}\n".format(map_read))
            for line in interval_report:
                report.write("\t".join(map(str, line))+"\n")
                
    def _write_bam(self, al_file, header, read_list):
        """
        Write all read mapped on intervals in a sorted bam file
        """
        # Generate a output name for the bam file
        outname = "./{}_overlap".format(file_basename(al_file))
        
        with pysam.Samfile("temp.bam", "wb", header=header) as bamfile:
            for read in read_list:
                bamfile.write(read)
        
        # That is a little dirty and unoptimized but it works
        pysam.sort("temp.bam", outname)
        os.remove("temp.bam")
        pysam.index (outname+".bam")


#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    main = Main()
    main()
