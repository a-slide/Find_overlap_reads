
# GLOBAL IMPORTS

import pysam # from pysam 0.8.1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Interval(object):
    """
    @class  Reference
    @brief  Object oriented class containing informations of genomic interval
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
    
    #~~~~~~~CLASS FIELDS~~~~~~~#

    Instances = [] # Class field used for instance tracking
    id_count = 1

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def next_id (self):
        cur_id = self.id_count
        self.id_count +=1
        return cur_id

    @ classmethod
    def countInstances (self):
        return len(self.Instances)

    @ classmethod
    def getInstances (self):
        return self.Instances

    @ classmethod
    def printInstances (self):
        for ref in self.Instances:
            print (repr(ref))

    @ classmethod
    def resetInstances (self):
        print ("Clearing Reference instances list")
        self.Instances = []
        self.id_count = 0

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, ref_name, start, end):
        """
        """
        # Store object variables
        self.id = self.next_id()
        self.ref_name = ref_name
        
        # Store start and end in crescent order
        if start <= end:
            self.start = int(start)
            self.end = int(end)
        else:
            self.end = int(start)
            self.start = int(end)

        # Define additional variables
        self.nread = 0
        self.read_list = []

        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __str__(self):
        return "{} [{}:{}] = {} overlapping read pairs".format(
            self.ref_name,
            self.start,
            self.end,
            self.nread)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __len__(self):
        return self.nread

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def is_overlapping (self, ref_name, start, end):
        
        # Reverse value order if negative order 
        if start > end:
            start, end = end, start
        
        return ref_name == self.ref_name and start <= self.start and end >= self.end
    
    def add_read (self, read):
        """
        Add a read to read_list and update the counter
        """
        self.read_list.append(read)
        self.nread+=1

    def sort_read (self):
        """
        sort read in read_list according to their leftmost position
        """
        self.read_list.sort(key = lambda x: x.pos)
        
    ## def WRITE BAM (self):
