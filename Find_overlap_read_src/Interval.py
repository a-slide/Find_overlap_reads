# GLOBAL IMPORTS
#import pysam # from pysam 0.8.1

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
        for inter in self.Instances:
            print (inter)

    @ classmethod
    def resetInstances (self):
        self.Instances = []
        self.id_count = 0
    
    @ classmethod
    def get_read (self):
        read_list =[]
        for i in self.Instances:
            read_list.extend(i.read_list)
        return read_list
    
    @ classmethod
    def get_report (self):
        report = [["ref_name", "start", "end", "name", "nread"]]
        report += [[i.ref_name, i.start, i.end, i.name, i.nread] for i in self.Instances]
        return report
    
    @ classmethod
    def resetReadCount (self):
        for inter in self.Instances:
            inter.nread=0
        
    @ classmethod
    def resetReadList (self):
        for inter in self.Instances:
            inter.read_list=[]

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, ref_name, start, end, name="-"):
        """
        @param ref_name Name of the reference sequence
        @param start Start coordinates of the interval (INT)
        @param end End coordinates of the interval (INT)
        @param name Facultative name of the interval
        """
        # Store object variables
        self.id = self.next_id()
        self.ref_name = ref_name
        self.name = name

        # Store start and end in crescent order
        if start <= end:
            self.start, self.end = start, end
        else:
            self.start, self.end = end, start

        # Define additional variables
        self.nread = 0
        self.read_list = []

        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __str__(self):
        return "{} [{}:{}] {} = {} reads found".format(
            self.ref_name,
            self.start,
            self.end,
            self.name,
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

