from collections import OrderedDict
import logging
import os
import sys

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph
from multiqc.base_module import ModuleNoSamplesFound

import yaml

log = logging.getLogger('damlabqc.strainline')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        super(MultiqcModule, self).__init__(
            name='Strainline',
            anchor='strainline',
            href='https://github.com/HaploKit/Strainline',
            info='is a pipeline for haplotype reconstruction of viral sequences.'
        )

        # Initialize data storage
        self._strainline_data = dict()
        
        self._collect_log_files()
        self._add_general_stats()
        self._add_sections()
    
    def _collect_log_files(self):
        """Collect data from log files found my multiqc."""
        
        for f in self.find_log_files('strainline'):
            self._parse_strainline_log(f)

        # Filter out ignored samples if we have any using multiqc public method
        self._strainline_data = self.ignore_samples(self._strainline_data)

        # No samples found
        if len(self._strainline_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self._strainline_data)} reports")
        
        # Write parsed report data to a file using multiqc public method
        self.write_data_file(self._strainline_data, 'multiqc_strainline')
    
    
    
    def _add_general_stats(self):
        """Add to the general stats table at the top of the report."""

        # Create general stats data from parsed strainline data
        general_stats_data = {}
        for s_name, metrics in self._strainline_data.items():
            general_stats_data[s_name] = {
                'haplotype_count': metrics['haplotype_count'],
                'haplotype_freqs_max': metrics['haplotype_freqs']['max']
            }

        # Make nice header columns
        headers = {
            'haplotype_count': {
                'title':'Haplotype Count',
                'description':'Number of haplotypes reconstructed',
                'scale':'RdYlGn',
                'format':'{:,.0f}',
                'placement':900.0  # Place towards left side of table
            },
            'haplotype_freqs_max': {
                'title':'Max Frequency',
                'description':'Maximum frequency of any haplotype',
                'scale':'RdYlGn',
                'min':0,
                'max':1,
                'format': '{:,.1%}',
                'placement':910.0
            }
        }
        
        # Add to general stats table using multiqc public method
        self.general_stats_addcols(general_stats_data, headers)
        

    def _parse_strainline_log(self, f):
        """Parse single strainline log file."""

        try:
            # Load YAML data
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            # Get sample name from the data
            s_name = data.get('sample_name', f['fn'])
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
            
            self._strainline_data[s_name] = metrics
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
        
    def _add_sections(self):
        """Add sections to MultiQC report."""
        
        self._add_haplotype_counts()
        self._add_haplotype_lengths()
        self._add_haplotype_freqs()
        
    def _add_haplotype_counts(self):
        """Add haplotype count barplot to report."""
        
        count_data = {}
        for s_name, d in self._strainline_data.items():
            count_data[s_name] = {'haplotype_count': d['haplotype_count']}

        self.add_section(
            name='Haplotype Counts',
            anchor='strainline-haplotype-counts',
            description='Number of haplotypes reconstructed per sample',
            plot=bargraph.plot(count_data,
                               pconfig = {
                                   'id': 'strainline_haplotype_counts',
                                   'title': 'Strainline: Haplotype Counts',
                                   'ylab': 'Number of Haplotypes'
                               }
            )
        )
    
    def _add_haplotype_lengths(self):
        """Add haplotype length distribution to report."""
        
        length_data = {
            s_name: {
                'Maximum': d['haplotype_lengths']['max'],
                'Mean': d['haplotype_lengths']['mean'],
                'Minimum': d['haplotype_lengths']['min']
            }
            for s_name, d in self._strainline_data.items()
        }
        
        self.add_section(
            name='Haplotype Lengths',
            anchor='strainline-haplotype-lengths',
            description='Distribution of haplotype lengths per sample',
            plot=bargraph.plot(length_data,
                               pconfig = {
                                   'id': 'strainline_haplotype_lengths',
                                   'title': 'Strainline: Haplotype Lengths',
                                   'ylab': 'Length (bp)'
                               }
            )
        )

    def _add_haplotype_freqs(self):
        """Add haplotype frequency distribution to report."""
        
        freq_data = {
            s_name: {
                'Maximum': d['haplotype_freqs']['max'],
                'Mean': d['haplotype_freqs']['mean'],
                'Minimum': d['haplotype_freqs']['min']
            }
            for s_name, d in self._strainline_data.items()
        }
        
        self.add_section(
            name='Haplotype Frequencies',
            anchor='strainline-haplotype-freqs',
            description='Distribution of haplotype frequencies per sample',
            plot=bargraph.plot(freq_data,
                               pconfig = {
                                   'id': 'strainline_haplotype_freqs',
                                   'title': 'Strainline: Haplotype Frequencies',
                                   'ylab': 'Frequency'
                               }
            )
        ) 