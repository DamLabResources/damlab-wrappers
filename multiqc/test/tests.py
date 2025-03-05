import pytest
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


def test_dorado():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "dorado-total_reads",
        "dorado-duplex_reads",
        "dorado-duplex_rate",
        "dorado-simplex_qscore",
        "dorado-duplex_qscore"
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data (5 reads: 2 duplex, 1 simplex with duplex, 1 simplex no duplex, 1 no tag)
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    assert float(sample1['dorado-total_reads']) == 5
    assert float(sample1['dorado-duplex_reads']) == 2
    assert float(sample1['dorado-duplex_rate']) == pytest.approx(0.4)  # 2/5
    assert float(sample1['dorado-simplex_qscore']) == pytest.approx(13.27, abs=0.1)  # Mean of 12.5, 14.2, 13.1
    assert float(sample1['dorado-duplex_qscore']) == pytest.approx(18.95, abs=0.1)  # Mean of 18.7, 19.2
    
    # Check sample2 data (6 reads: 2 duplex, 2 simplex with duplex, 1 simplex no duplex, 1 no tag)
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['dorado-total_reads']) == 6
    assert float(sample2['dorado-duplex_reads']) == 2
    assert float(sample2['dorado-duplex_rate']) == pytest.approx(0.333, abs=0.01)  # 2/6
    
    # Check sample3 data (7 reads: 3 duplex, 2 simplex with duplex, 1 simplex no duplex, 1 no tag)
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['dorado-total_reads']) == 7
    assert float(sample3['dorado-duplex_reads']) == 3
    assert float(sample3['dorado-duplex_rate']) == pytest.approx(0.429, abs=0.01)  # 3/7


def test_generic():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that the wanted column is present
    assert 'generic-gel_band' in data[0]
    
    # Check that the data is correct for each sample
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    
    assert sample1['generic-gel_band'] == 'yes'
    assert sample2['generic-gel_band'] == 'yes'
    assert sample3['generic-gel_band'] == 'no'


def test_primercheck():
    data = get_general_stats_data()

    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "primercheck-total_sequences",
        "primercheck-primer_primer1",
        "primercheck-primer_primer2",
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data (4/5 reads with primer1, 1/5 with primer2/3)
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    assert float(sample1['primercheck-total_sequences']) == 5
    assert float(sample1['primercheck-primer_primer1']) == pytest.approx(0.8)  # 4/5 reads
    assert float(sample1['primercheck-primer_primer2']) == pytest.approx(0.2)  # 1/5 reads
    
    # Check sample2 data (2/5 with primer1, 3/5 with primer2/3)
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['primercheck-total_sequences']) == 5
    assert float(sample2['primercheck-primer_primer1']) == pytest.approx(0.4)  # 2/5 reads
    assert float(sample2['primercheck-primer_primer2']) == pytest.approx(0.6)  # 3/5 reads
    
    # Check sample3 data (all reads with primer2/3)
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['primercheck-total_sequences']) == 4
    assert float(sample3['primercheck-primer_primer1']) == pytest.approx(0.0)  # No reads
    assert float(sample3['primercheck-primer_primer2']) == pytest.approx(1.0)  # All reads

