__author__ = "Will Dampier"
__copyright__ = "Copyright 2022, Will Dampier"
__email__ = "wnd22@drexel.edu"
__license__ = "Kept"

import pandas as pd
import logomaker as lm

from Bio import motifs


def get_motif(path, name):
    """Gets the motif from a jaspar formatted file that matched the name.
    If the name is None, returns the first"""

    with open(path) as handle:
        for motif in motifs.parse(handle, 'jaspar'):
            if name is None:
                return motif
            elif name == motif.matrix_id:
                return motif

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake

# inputs

jaspar_path = snakemake.input.get("jaspar")

# output:

logo_path = snakemake.output.get("logo")

# params:

dpi = snakemake.params.get("dpi")
name_filter = snakemake.params.get("name_filter")
target_sequence = snakemake.params.get("target_sequence")

motif = get_motif(jaspar_path, name_filter)

if motif is not None:
    pwm_df = pd.DataFrame(motif.pwm)

    inf_df = lm.transform_matrix(pwm_df,
                                 from_type="probability",
                                 to_type="information")
else:
    inf_df = lm.sequence_to_matrix('N'*len(target_sequence),
                                   is_iupac=True,
                                   to_type='information')
    

logo = lm.Logo(inf_df)
logo.style_spines(spines=["top", "right"], visible=False)
logo.highlight_position_range(pmin=22, pmax=27, color="lightcyan")
logo.ax.set_ylabel("Information")
logo.ax.set_title(name_filter)

if target_sequence is not None:

    logo.ax.set_xticks(range(len(target_sequence)))
    logo.ax.set_xticklabels(list(target_sequence))

logo.fig.savefig(logo_path, dpi=dpi)
