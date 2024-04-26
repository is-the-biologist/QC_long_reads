import subprocess
from qc_module import library_statistics
import numpy as np
import pandas as pd
import os


#Unsophisticated wrapper to run minimal QC methods with only a txt with filepaths to BAMs as input
if __name__ == '__main__':

    npnames = []
    with open('filepaths.txt', 'r') as myFiles:
        for line in myFiles:
            bam_file = line.strip("\n")

            if os.path.isdir("readlengths"): #make dir to store numpy files to make histograms later
                pass
            else:
                os.mkdir("readlengths")

            output = "readlengths/"+bam_file.split('/')[-1].replace("bam","npy")
            ps = subprocess.Popen(["samtools", "view", "-h", f"{bam_file}"], stdout=subprocess.PIPE)
            ps2 = subprocess.Popen(('cut', '-f10'), stdin=ps.stdout, stdout=subprocess.PIPE)
            ps3 = subprocess.check_output(("python", "generateReadLengths.py", f"{output}"), stdin=ps2.stdout)
            ps.wait()
            npnames.append(output)
        myFiles.close()

    stat_matrix = []
    libstats = library_statistics()

    for input_library in npnames:
        RL_np = np.load(input_library)
        stats = libstats.generateStatsTable(RL_np)
        libstats.histogram(RL_array=RL_np, snk_out=input_library)
        stat_matrix.append(stats)

    stat_table = pd.DataFrame(data=stat_matrix,
                              columns=['Total_bp', 'Mean', "SD", "1-ile", '10-ile', '25-ile', '50-ile', '75-ile', '90-ile', '99-ile'],
                              index=[fname.split('/')[-1].replace(".npy", "") for fname in npnames])
    stat_table.to_csv("library_stats.csv", index=True)