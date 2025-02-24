from collections import OrderedDict
import logging
import os
import yaml
import numpy as np

from multiqc import config
from multiqc.plots import bargraph, linegraph, heatmap
from multiqc.base_module import BaseMultiqcModule
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
        """Add read length distribution plot as a heatmap showing percentage of reads at each length."""
        # Prepare data structures
        bin_centers = []
        samples = []
        
        # First pass to get bin centers and samples
        for s_name, d in self._hivmetrics_data.items():
            if d['length_histogram'] is not None:
                hist = d['length_histogram']
                if not bin_centers:  # Only need to do this once
                    bin_centers = [(a + b)/2 for a, b in zip(hist['bin_edges'][:-1], hist['bin_edges'][1:])]
                samples.append(s_name)
        
        # Initialize normalized matrix with zeros (transposed orientation)
        num_bins = len(bin_centers)
        num_samples = len(samples)
        normalized_data = [[0] * num_samples for _ in range(num_bins)]
        
        # Fill matrix (transposed)
        for sample_idx, s_name in enumerate(samples):
            hist = self._hivmetrics_data[s_name]['length_histogram']
            raw_counts = hist['counts']
            total_reads = sum(raw_counts)
            
            for length_idx, count in enumerate(raw_counts):
                # Normalized counts (percentage of reads for this sample)
                if total_reads > 0:
                    normalized_data[length_idx][sample_idx] = count/total_reads
        
        # Reverse the order of rows so longer fragments are at the top
        normalized_data.reverse()
        bin_centers.reverse()

        if normalized_data:
            # Format length labels for y-axis
            y_labels = [f"{int(x):,}" for x in bin_centers]
            
            self.add_section(
                name='Read Length Distribution',
                anchor='hivmetrics-length-dist',
                description='Distribution of mapped read lengths. Each column represents a sample, and the color intensity indicates the percentage of reads at each length. Longer fragments appear near the top, similar to an agarose gel.',
                plot=heatmap.plot(
                    normalized_data,
                    xcats=samples,
                    ycats=y_labels,
                    pconfig={
                        'id': 'hivmetrics_length_dist',
                        'title': 'HIVmetrics: Read Length Distribution',
                        'xlab': 'Sample',
                        'ylab': 'Read Length (bp)',
                        'zlab': 'Percentage of Reads',
                        'min': 0,
                        'max': 1,
                        'square': False,
                        'colstops': [
                            [0, '#FFFFFF'],
                            [0.001, '#FFF7EC'],  # Very light color at 0.1%
                            [0.01, '#FEE8C8'],   # Light color at 1%
                            [0.1, '#FDD49E'],    # Medium color at 10%
                            [0.5, '#FC8D59'],    # Darker color at 50%
                            [1, '#B30000']       # Deep red at 100%
                        ],
                        'height': 500,
                        'ycats_samples': False,
                        'xcats_samples': True
                    }
                )
            ) 