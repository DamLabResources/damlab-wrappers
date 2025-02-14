import pytest # type: ignore
import os
import csv


def test_multiqc_report_exists():
    assert os.path.exists("build/multiqc/multiqc_report.html")


def test_multiqc_report_general_stats_exists():
    assert os.path.exists("build/multiqc/multiqc_report_data/multiqc_general_stats.txt")

def get_general_stats_data():
    with open("build/multiqc/multiqc_report_data/multiqc_general_stats.txt", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def test_strainline():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that the wanted columns are present
    wanted_cols = ["strainline-haplotype_count", "strainline-haplotype_freqs_max"]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check that the data is correct
    assert [row['strainline-haplotype_count'] for row in data] == ['3', '1', '2']
    assert [row['strainline-haplotype_freqs_max'] for row in data] == ['0.75', '1.0', '0.5']

