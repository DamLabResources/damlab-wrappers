from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule # type: ignore
from multiqc.plots import linegraph, bargraph # type: ignore
from multiqc.plots.table_object import ColumnMeta # type: ignore
import logging
import os
import yaml # type: ignore
from collections import OrderedDict
import sys

log = logging.getLogger('multiqc.modules.strainline')
log.critical(f"Got log somehow")

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        

        super(MultiqcModule, self).__init__(
            name='Strainline',
            anchor='strainline',
            href='https://github.com/HaploKit/Strainline',
            info='is a pipeline for haplotype reconstruction of viral sequences.'
        )

        # Add to module config
        headers = {
            'haplotype_count': ColumnMeta(
                title='Haplotype Count',
                description='Number of haplotypes reconstructed',
                scale='RdYlGn',
                format='{:,.0f}',
                placement=900.0  # Place towards left side of table
            ),
            'haplotype_freqs_max': ColumnMeta(
                title='Max Frequency',
                description='Maximum frequency of any haplotype',
                scale='RdYlGn',
                min=0,
                max=1,
                format='{:,.3f}',
                placement=910.0
            )
        }
        
        config.module_order = ["strainline"] + config.module_order
        config.general_stats_headers['strainline'] = headers

        # Find all strainline output files
        self.strainline_data = dict()
        
        log.critical(f"Checking for strainline log files")
        # Look for strainline log files
        for f in self.find_log_files('strainline'):
            log.critical(f"Found strainline log file: {f['fn']}")
            
            self.parse_strainline_logs(f)
            
        # Filter samples if we have any
        self.strainline_data = self.ignore_samples(self.strainline_data)
        
        # No samples found
        if len(self.strainline_data) == 0:
            raise UserWarning
            
        log.info(f"Found {len(self.strainline_data)} reports")
        
        # Write parsed report data to a file
        self.write_data_file(self.strainline_data, 'multiqc_strainline')

        # Add to general stats table
        general_stats_data = {}
        for s_name, metrics in self.strainline_data.items():
            general_stats_data[s_name] = {
                'haplotype_count': metrics['haplotype_count'],
                'haplotype_freqs_max': metrics['haplotype_freqs']['max']
            }
        self.general_stats_addcols(general_stats_data, headers)
        
        # Add sections to report
        self.add_sections()
        
    def parse_strainline_logs(self, f):
        """Parse strainline log files."""
        try:
            # Load YAML data
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            # Get sample name from the data
            s_name = data.get('sample_id', f['s_name'])
            s_name = self.clean_s_name(s_name, f)
            
            # Parse metrics
            metrics = {
                'haplotype_count': data.get('haplotype_count', 0),
                'haplotype_lengths': {
                    'max': data.get('haplotype_max_length', 0),
                    'min': data.get('haplotype_min_length', 0),
                    'mean': data.get('haplotype_mean_length', 0)
                },
                'haplotype_freqs': {
                    'max': data.get('haplotype_max_freq', 0),
                    'min': data.get('haplotype_min_freq', 0),
                    'mean': data.get('haplotype_mean_freq', 0)
                }
            }
            
            self.strainline_data[s_name] = metrics
            
        except Exception as e:
            log.debug(f"Error parsing file {f['fn']}: {e}")
            return
        
    def add_sections(self):
        """Add sections to MultiQC report."""
        # Add haplotype count barplot
        self.add_section(
            name='Haplotype Counts',
            anchor='strainline-haplotype-counts',
            description='Number of haplotypes reconstructed per sample',
            plot=bargraph.plot({
                'id': 'strainline_haplotype_counts',
                'data': {s_name: {'haplotype_count': d['haplotype_count']}
                        for s_name, d in self.strainline_data.items()},
                'title': 'Strainline: Haplotype Counts',
                'ylab': 'Number of Haplotypes'
            })
        )
        
        # Add haplotype length distribution
        length_data = {
            s_name: {
                'Maximum': d['haplotype_lengths']['max'],
                'Mean': d['haplotype_lengths']['mean'],
                'Minimum': d['haplotype_lengths']['min']
            }
            for s_name, d in self.strainline_data.items()
        }
        
        self.add_section(
            name='Haplotype Lengths',
            anchor='strainline-haplotype-lengths',
            description='Distribution of haplotype lengths per sample',
            plot=bargraph.plot({
                'id': 'strainline_haplotype_lengths',
                'data': length_data,
                'title': 'Strainline: Haplotype Lengths',
                'ylab': 'Length (bp)'
            })
        )
        
        # Add haplotype frequency distribution
        freq_data = {
            s_name: {
                'Maximum': d['haplotype_freqs']['max'],
                'Mean': d['haplotype_freqs']['mean'],
                'Minimum': d['haplotype_freqs']['min']
            }
            for s_name, d in self.strainline_data.items()
        }
        
        self.add_section(
            name='Haplotype Frequencies',
            anchor='strainline-haplotype-freqs',
            description='Distribution of haplotype frequencies per sample',
            plot=bargraph.plot({
                'id': 'strainline_haplotype_freqs',
                'data': freq_data,
                'title': 'Strainline: Haplotype Frequencies',
                'ylab': 'Frequency'
            })
        ) 