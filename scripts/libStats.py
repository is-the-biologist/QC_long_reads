from qc_module import library_statistics
import numpy as np
import pandas as pd

if __name__ == '__main__':
    stat_matrix = []
    libstats = library_statistics()
    for input_library in snakemake.input:
        RL_np = np.load(input_library)
        stats = libstats.generateStatsTable(RL_np)
        libstats.histogram(RL_array=RL_np, snk_out=input_library)
        stat_matrix.append(stats)

    stat_table = pd.DataFrame(data=stat_matrix,
                              columns=['Total bp', 'Mean', "SD", "1-ile", '10-ile', '25-ile', '50-ile', '75-ile', '90-ile', '99-ile'],
                              index=[fname.split('/')[-1] for fname in snakemake.input])
    stat_table.to_csv(snakemake.output[0], index=True)
