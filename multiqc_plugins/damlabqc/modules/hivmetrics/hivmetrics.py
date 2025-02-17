from collections import OrderedDict
import logging
import os
import yaml
import numpy as np

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.hivmetrics')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='HIVmetrics',
            anchor='hivmetrics',
            href='',
            info='analyzes HIV sequencing metrics.'
        )

        # Initialize data storage
        self._hivmetrics_data = dict()
        
        self._collect_log_files()
        self._add_general_stats()
        self._add_sections()
    
    def _collect_log_files(self):
        """Collect data from log files found by multiqc."""
        for f in self.find_log_files('hivmetrics'):
            self._parse_hivmetrics_log(f)

        # Filter out ignored samples
        self._hivmetrics_data = self.ignore_samples(self._hivmetrics_data)

        # No samples found
        if len(self._hivmetrics_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self._hivmetrics_data)} reports")
        
        # Write parsed report data
        self.write_data_file(self._hivmetrics_data, 'multiqc_hivmetrics')
    
    def _parse_hivmetrics_log(self, f):
        """Parse single hivmetrics log file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            
            self._hivmetrics_data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add to the general stats table at the top of the report."""
        general_stats_data = {}
        
        for s_name, data in self._hivmetrics_data.items():
            general_stats_data[s_name] = {
                'total_reads': data['total_reads'],
                'mapped_reads': data['mapped_reads'],
                'mapping_rate': data['mapping_rate'],
                'supplementary_reads': data['supplementary_reads'],
                'mean_length': data['mean_mapped_length']
            }

        headers = {
            'total_reads': {
                'title': 'Total Reads',
                'description': 'Total number of reads',
                'scale': 'GnBu',
                'format': '{:,.0f}',
                'placement': 1000
            },
            'mapped_reads': {
                'title': 'Mapped',
                'description': 'Number of mapped reads',
                'scale': 'YlGn',
                'format': '{:,.0f}',
                'placement': 1100
            },
            'mapping_rate': {
                'title': 'Mapping Rate',
                'description': 'Proportion of reads that mapped',
                'max': 1,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.1%}',
                'placement': 1200
            },
            'supplementary_reads': {
                'title': 'Supplementary',
                'description': 'Number of supplementary alignments',
                'scale': 'PuRd',
                'format': '{:,.0f}',
                'placement': 1300
            },
            'mean_length': {
                'title': 'Mean Length',
                'description': 'Mean length of mapped reads',
                'scale': 'BuPu',
                'format': '{:,.0f}',
                'placement': 1400
            }
        }
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_sections(self):
        """Add sections to MultiQC report."""
        self._add_mapping_stats()
        self._add_length_distribution()
    
    def _add_mapping_stats(self):
        """Add mapping statistics plot."""
        mapping_data = {}
        
        for s_name, d in self._hivmetrics_data.items():
            mapping_data[s_name] = {
                'Mapped': d['mapped_reads'],
                'Unmapped': d['total_reads'] - d['mapped_reads'],
                'Supplementary': d['supplementary_reads']
            }

        self.add_section(
            name='Mapping Statistics',
            anchor='hivmetrics-mapping-stats',
            description='Distribution of read mapping results',
            plot=bargraph.plot(mapping_data,
                             pconfig={
                                 'id': 'hivmetrics_mapping_stats',
                                 'title': 'HIVmetrics: Mapping Statistics',
                                 'ylab': 'Number of Reads',
                                 'stacking': 'normal'
                             })
        )
    
    def _add_length_distribution(self):
        """Add read length distribution plot."""
        plot_data = {}
        
        for s_name, d in self._hivmetrics_data.items():
            if d['length_histogram'] is not None:
                hist = d['length_histogram']
                # Create x-coordinates as bin centers
                x = [(a + b)/2 for a, b in zip(hist['bin_edges'][:-1], hist['bin_edges'][1:])]
                plot_data[s_name] = {x_coord: y_val for x_coord, y_val in zip(x, hist['counts'])}

        if plot_data:
            self.add_section(
                name='Read Length Distribution',
                anchor='hivmetrics-length-dist',
                description='Distribution of mapped read lengths',
                plot=linegraph.plot(plot_data,
                                  pconfig={
                                      'id': 'hivmetrics_length_dist',
                                      'title': 'HIVmetrics: Read Length Distribution',
                                      'xlab': 'Read Length (bp)',
                                      'ylab': 'Number of Reads',
                                      'smooth_points': 200
                                  })
            ) 