#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@package    
@brief      
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

    
# IMPORTS
try:

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()

#~~~~~~~MAIN FUNCTION~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
    """
    
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self):

        if len(argv) != 2:
            print ("Please provide the path to the Configuration file as an unique argument\n")
            exit()


    def __repr__(self):
        msg = "MAIN CLASS\n"
        msg+= "\tParameters list\n"
        for i, j in self.__dict__.items():
            msg+="\t{}\t{}\n".format(i, j)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        Launch the complete pipeline of analyse:

        * Reference importation/parsing
        * Facultative step of reference masking to remove homologies between reference sequences
        * Facultative step of Fastq quality Filtering/ adapter trimming
        * Facultative step of reference indexing for bwa from merged references
        * Short read alignment with bwa mem
        * Spliting of sam to attribute reads to each original references (or unmmapped)
        * Output per reference bam, sam, bedgraph, bed, covgraph, variant call
        * Output distribution table and graph
        """
        stime = time()



    ##~~~~~~~PRIVATE METHODS~~~~~~~#

    def _optparser(self, program_name, program_version):
        """
        Parse command line arguments with optparse and verify the file validity
        @param program_name Name of the program
        @param program_version Version of the program
        @return A dictionnary containing:
        * Host and virus fasta file paths
        * Configuration file path
        * Sequencing mode : pair or single
        * Basename for output files
        """
        # Usage and version strings
        usage_string = "%prog -H Host_genome.fa[.gz] -V Viral_genome.fa[.gz] -C Conf_file.txt [-o Output_prefix] [-p |-s]"
        version_string = program_name + program_version
        optparser = optparse.OptionParser(usage = usage_string, version = version_string)

        # Define optparser options
        hstr = "Path of the fasta file containing the host genome sequence (can be gziped)"
        optparser.add_option( '-H', '--host_genome', dest="hg", help=hstr)
        hstr = "Path of the fasta file containing the viral genome sequence (can be gziped)"
        optparser.add_option( '-V', '--virus_genome', dest="vg", help=hstr)
        hstr = "Path of the configuration text file"
        optparser.add_option( '-C', '--conf_file', dest="conf", help=hstr)
        hstr = "Facultative option to indicate the name of the output prefix (default = out)"
        optparser.add_option( '-o', '--output', default="out", dest="output", help=hstr)
        htr = "Single end mode overwriten if -p option is also indicated"
        optparser.add_option( '-s', '--single', dest="single", action='store_true', help=hstr)
        hstr = "Pair end mode incompatible with -s option (default mode)"
        optparser.add_option( '-p', '--pair', dest="pair", action='store_true', help=hstr)

        # Parse arg and return a dictionnary_like object of options
        options, args = optparser.parse_args()

        # Validate option and generate a dictionnary
        arg_dict = {'host_genome' : self._check_file (options.hg, "host_genome"),
                    'virus_genome' : self._check_file (options.vg, "virus_genome "),
                    'conf_file' : self._check_file (options.conf, "conf_file"),
                    'basename' : options.output,
                    'pair' : self._check_mode (options.single, options.pair)}

        return arg_dict



	def _check_file (self, path, descr):
        """
        Try to open the file at a given path. In case of impossibility an IsisConfException is raised
        @param path     Path of the file to verify
        @param descr    name of the option associated with the file
        @return         Valid path
        """
        if not path:
            raise IsisConfException ("{} is a mandatory command line argument".format(descr))

        try:
            handle = open(path, "r")
            handle.close()
            print ("\t\tValid file for {}".format (descr))
            return path

        except IOError:
            raise IsisConfException ("Error : " + path + " can't be read\n\
            Please enter a valid path")


    def _sam_spliter (self):
        """
        """
        with pysam.Samfile(self.sam, "r" ) as samfile:
            self.bam_header = samfile.header

            # Give the header of the sam file to all Reference.Instances to respect the same order
            # references in sorted bam files
            Reference.set_global("bam_header", self.bam_header)

            # Create a dict to collect unmapped and low quality reads
            Secondary = Sequence (name = 'Secondary', length = 0)
            Unmapped = Sequence (name = 'Unmapped', length = 0)
            LowMapq = Sequence (name = 'LowMapq', length = 0)
            self.garbage_read = [Secondary, Unmapped, LowMapq]

            for read in samfile:
                # Always remove secondary alignments
                if read.is_secondary:
                    Secondary.add_read(read)
                # Filter Unmapped reads
                elif read.tid == -1:
                    Unmapped.add_read(read)
                # Filter Low MAPQ reads
                elif read.mapq < self.min_mapq:
                    LowMapq.add_read(read)
                # Finally if all is fine attribute the read to a Reference
                else:
                    Reference.addRead(samfile.getrname(read.tid), read)
        
        # Removing the original sam file which is no longer needed
        remove(self.sam)
        self.sam = None






#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    main = Main()
    main()	
