from qc_module import library_statistics
import sys
import numpy as np
if __name__ == '__main__':
    snkmk_out = sys.argv[1]

    libcall = library_statistics()
    read_length_array = np.fromiter(libcall.read_in_bam_stream(), dtype="uint32")

    np.save(snkmk_out, read_length_array)