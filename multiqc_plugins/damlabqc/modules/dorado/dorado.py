from collections import OrderedDict
import logging
import os
import yaml

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.dorado')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Dorado',
            anchor='dorado',
            href='https://github.com/nanoporetech/dorado',
            info='is a basecaller for Oxford Nanopore sequencing data.'
        )

        # Initialize data storage
        self._dorado_data = dict()
        
        self._collect_log_files()
        self._add_general_stats()
        self._add_sections()
    
    def _collect_log_files(self):
        """Collect data from log files found by multiqc."""
        for f in self.find_log_files('dorado'):
            self._parse_dorado_log(f)

        # Filter out ignored samples
        self._dorado_data = self.ignore_samples(self._dorado_data)

        # No samples found
        if len(self._dorado_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self._dorado_data)} reports")
        
        # Write parsed report data
        self.write_data_file(self._dorado_data, 'multiqc_dorado')
    
    def _parse_dorado_log(self, f):
        """Parse single dorado log file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            
            self._dorado_data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add to the general stats table at the top of the report."""
        general_stats_data = {}
        
        for s_name, data in self._dorado_data.items():
            general_stats_data[s_name] = {
                'total_reads': data['simplex']['total_reads'] + data['duplex']['total_reads'],
                'duplex_reads': data['duplex']['total_reads'],
                'duplex_rate': data['duplex']['total_reads'] / (data['simplex']['total_reads'] + data['duplex']['total_reads']) if (data['simplex']['total_reads'] + data['duplex']['total_reads']) > 0 else 0,
                'simplex_qscore': data['simplex']['mean_qscore'],
                'duplex_qscore': data['duplex']['mean_qscore']
            }

        headers = {
            'total_reads': {
                'title': 'Total Reads',
                'description': 'Total number of reads (simplex + duplex)',
                'scale': 'GnBu',
                'format': '{:,.0f}',
                'placement': 1000
            },
            'duplex_reads': {
                'title': 'Duplex Reads',
                'description': 'Number of duplex reads',
                'scale': 'PuRd',
                'format': '{:,.0f}',
                'placement': 1100
            },
            'duplex_rate': {
                'title': 'Duplex Rate',
                'description': 'Proportion of reads that are duplex',
                'max': 1,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.1%}',
                'placement': 1200
            },
            'simplex_qscore': {
                'title': 'Simplex Q',
                'description': 'Mean quality score for simplex reads',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.1f}',
                'placement': 1300
            },
            'duplex_qscore': {
                'title': 'Duplex Q',
                'description': 'Mean quality score for duplex reads',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.1f}',
                'placement': 1400
            }
        }
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_sections(self):
        """Add sections to MultiQC report."""
        self._add_read_counts()
        self._add_read_lengths()
        #self._add_pass_rates()
    
    def _add_read_counts(self):
        """Add read count distribution plots."""
        reads_data = {}
        bases_data = {}
        
        for s_name, d in self._dorado_data.items():
            reads_data[s_name] = {
                'Simplex (No Duplex)': d['simplex']['no_duplex'],
                'Simplex (Has Duplex)': d['simplex']['has_duplex'],
                'Duplex': d['duplex']['total_reads']
            }
            bases_data[s_name] = {
                'Simplex': d['simplex']['total_bases'],
                'Duplex': d['duplex']['total_bases']
            }

        self.add_section(
            name='Read Counts',
            anchor='dorado-read-counts',
            description='Distribution of read types and bases',
            plot=bargraph.plot([reads_data, bases_data],
                             pconfig={
                                 'id': 'dorado_read_counts',
                                 'title': 'Dorado: Read & Base Counts',
                                 'ylab': 'Count',
                                 'data_labels': ['Reads', 'Bases'],
                                 'stacking': 'normal'
                             })
        )
    
    def _add_read_lengths(self):
        """Add read length metrics plot."""
        length_data = {}
        
        for s_name, d in self._dorado_data.items():
            length_data[s_name] = {
                'Simplex Mean': d['simplex']['mean_read_length'],
                'Simplex N50': d['simplex']['read_length_n50'],
                'Duplex Mean': d['duplex']['mean_read_length'],
                'Duplex N50': d['duplex']['read_length_n50']
            }

        self.add_section(
            name='Read Lengths',
            anchor='dorado-read-lengths',
            description='Read length statistics by read type',
            plot=bargraph.plot(length_data,
                             pconfig={
                                 'id': 'dorado_read_lengths',
                                 'title': 'Dorado: Read Lengths',
                                 'ylab': 'Length (bp)'
                             })
        )
    
    def _add_pass_rates(self):
        """Add pass rate comparison plot."""
        pass_data = {}
        
        for s_name, d in self._dorado_data.items():
            pass_data[s_name] = {
                'Simplex': d['simplex']['pass_rate'],
                'Duplex': d['duplex']['pass_rate']
            }

        self.add_section(
            name='Pass Rates',
            anchor='dorado-pass-rates',
            description='Pass rates by read type',
            plot=bargraph.plot(pass_data,
                             pconfig={
                                 'id': 'dorado_pass_rates',
                                 'title': 'Dorado: Pass Rates',
                                 'ylab': 'Pass Rate',
                                 'max': 1,
                                 'min': 0,
                                 'format': '{:.1%}'
                             })
        )