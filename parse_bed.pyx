from collections import defaultdict
import pandas as pd
import os
from sequencing_tools.io_tools import xopen

cpdef int get_length(str fragment):
    cdef:
        long start, end

    fields = fragment.split('\t')
    start, end = map(long, [fields[1], fields[2]])
    return abs(end - start)


def parse_bed(bed_file):
    cdef:
        str fragment
        int fragment_size

    samplename = os.path.basename(bed_file)
    print 'Extracting %s' %samplename
    isize_dict = defaultdict(int)

    with xopen(bed_file,'r') as bed:
        for fragment in bed:
            fragment_size = get_length(fragment)
            isize_dict[fragment_size] += 1

    df = pd.DataFrame({'isize': list(isize_dict.keys()),
                        'counts':list(isize_dict.values())})
    df['samplename'] = samplename

    temp_file = bed_file + '_temp'
    df.to_csv(temp_file, index=False)
    return temp_file
