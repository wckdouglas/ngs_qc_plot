from collections import defaultdict
import pandas as pd
import os

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

    samplename = os.path.basename(bed_file).split('.')[0]
    print 'Extracting %s' %samplename
    isize_dict = defaultdict(int)

    with open(bed_file,'r') as bed:
        for fragment in bed:
            fragment_size = get_length(fragment)
            isize_dict[fragment_size] += 1

    df = pd.DataFrame({'isize': isize_dict.keys(),
                        'counts':isize_dict.values()})
    df['samplename'] = samplename
    return df
