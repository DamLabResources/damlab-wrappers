from collections import OrderedDict, defaultdict
import logging
import yaml
import numpy as np

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.barcode')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Barcode',
            anchor='barcode',
            href='',
            info='analyzes barcode extraction and correction metrics.'
        )

        # Initialize data storage
        self.extract_data = dict()
        self.correct_data = dict()
        
        # Find and parse input files
        for f in self.find_log_files('barcode/extract'):
            self._parse_extract_yaml(f)
        
        for f in self.find_log_files('barcode/correct'):
            self._parse_correct_yaml(f)
            
        # Filter out ignored samples
        self.extract_data = self.ignore_samples(self.extract_data)
        self.correct_data = self.ignore_samples(self.correct_data)

        # No samples found
        if len(self.extract_data) == 0 and len(self.correct_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self.extract_data)} extraction samples")
        log.info(f"Found {len(self.correct_data)} correction samples")
        
        # Write parsed data to file
        self.write_data_file(self.extract_data, 'multiqc_barcode_extract')
        self.write_data_file(self.correct_data, 'multiqc_barcode_correct')
        
        # Add sections
        self._add_general_stats()
        self._add_extraction_plot()
        self._add_correction_plot()
        #self._add_count_distribution()
    
    def _parse_extract_yaml(self, f):
        """Parse extraction metrics YAML file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                return
                
            # Get sample name from data or filename
            s_name = data.get('sample_name', f['s_name'].split('.')[0])
            s_name = self.clean_s_name(s_name, f)
            
            self.extract_data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing extraction file {f['fn']}: {e}")
    
    def _parse_correct_yaml(self, f):
        """Parse correction metrics YAML file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                return
                
            # Combine sample_name and barcode_name for unique identifier
            s_name = data.get('sample_name', f['s_name'].split('.')[0])
            b_name = data.get('barcode_name', 'default')
            combined_name = f"{s_name}_{b_name}"
            combined_name = self.clean_s_name(combined_name, f)
            
            # Store original sample name for grouping
            data['_sample_name'] = s_name
            self.correct_data[combined_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing correction file {f['fn']}: {e}")
    
    def _add_general_stats(self):
        """Add overview stats to the general statistics table."""
        general_stats = {}
        
        # Process extraction data
        for s_name, d in self.extract_data.items():
            total = d.get('total_reads', 0)
            extraction_counts = d.get('extraction_counts', {})
            pairwise_counts = d.get('pairwise_counts', {})
            
            general_stats[s_name] = {
                'total_reads_extract': total,
                'barcode_rate': extraction_counts.get('barcode', 0) / total if total > 0 else 0,
                'all_extracted_rate': pairwise_counts.get('all_three', pairwise_counts.get('all', 0)) / total if total > 0 else 0
            }
        
        # Process correction data
        for s_name, d in self.correct_data.items():
            counts = d.get('barcode_counts', {})
            total_reads = counts.get('total_reads', 0)
            unique_before = counts.get('unique_barcodes_before', 0)
            correction_details = d.get('correction_details', {})
            barcode_name = d.get('barcode_name', 'default')
            
            # Get original sample name without barcode suffix
            orig_name = d.get('_sample_name', s_name.rsplit('_', 1)[0])
            
            # Create column names with barcode prefix
            col_prefix = f"{barcode_name}_" if barcode_name != 'default' else ''
            
            general_stats[orig_name] = general_stats.get(orig_name, {})
            general_stats[orig_name].update({
                f'{col_prefix}total_reads_correct': total_reads,
                f'{col_prefix}unique_before': unique_before,
                f'{col_prefix}unique_after': counts.get('unique_barcodes_after', 0),
                f'{col_prefix}correction_rate': correction_details.get('barcodes_corrected', 0) / unique_before if unique_before > 0 else 0
            })

        # Add columns to general stats table
        headers = OrderedDict()
        
        if self.extract_data:
            headers.update({
                'total_reads_extract': {
                    'title': 'Extract Reads',
                    'description': 'Total reads processed in extraction',
                    'format': '{:,.0f}',
                    'scale': 'Blues'
                },
                'barcode_rate': {
                    'title': 'Barcode Rate',
                    'description': 'Proportion of reads with extracted barcodes',
                    'max': 1,
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:.1%}'
                },
                'all_extracted_rate': {
                    'title': 'Full Extract Rate',
                    'description': 'Proportion of reads with all elements extracted',
                    'max': 1,
                    'min': 0,
                    'scale': 'RdYlGn',
                    'format': '{:.1%}'
                }
            })
        
        # Collect all unique barcode names
        barcode_names = set()
        for d in self.correct_data.values():
            barcode_names.add(d.get('barcode_name', 'default'))
            
        # Add headers for each barcode type
        if self.correct_data:
            for barcode_name in barcode_names:
                prefix = f"{barcode_name}_" if barcode_name != 'default' else ''
                headers.update({
                    f'{prefix}total_reads_correct': {
                        'title': f'{barcode_name} Reads' if barcode_name != 'default' else 'Correct Reads',
                        'description': f'Total reads processed in correction for {barcode_name}',
                        'format': '{:,.0f}',
                        'scale': 'Blues'
                    },
                    f'{prefix}unique_before': {
                        'title': f'{barcode_name} Initial' if barcode_name != 'default' else 'Initial Barcodes',
                        'description': f'Number of unique {barcode_name} barcodes before correction',
                        'format': '{:,.0f}',
                        'scale': 'Blues'
                    },
                    f'{prefix}unique_after': {
                        'title': f'{barcode_name} Final' if barcode_name != 'default' else 'Final Barcodes',
                        'description': f'Number of unique {barcode_name} barcodes after correction',
                        'format': '{:,.0f}',
                        'scale': 'Blues'
                    },
                    f'{prefix}correction_rate': {
                        'title': f'{barcode_name} Corr Rate' if barcode_name != 'default' else 'Correction Rate',
                        'description': f'Proportion of {barcode_name} barcodes that were corrected',
                        'max': 1,
                        'min': 0,
                        'scale': 'RdYlGn',
                        'format': '{:.1%}'
                    }
                })
        
        self.general_stats_addcols(general_stats, headers)
    
    def _add_extraction_plot(self):
        """Create barplot showing extraction success rates."""
        if not self.extract_data:
            return
            
        plot_data = {}
        for s_name, d in self.extract_data.items():
            extraction_counts = d.get('extraction_counts', {})
            pairwise_counts = d.get('pairwise_counts', {})
            
            plot_data[s_name] = {
                'Barcode': extraction_counts.get('barcode', 0),
                'Left UMI': extraction_counts.get('left_umi', 0),
                'Right UMI': extraction_counts.get('right_umi', 0),
                'All Three': pairwise_counts.get('all', 0)
            }
            
        self.add_section(
            name='Extraction Success',
            anchor='barcode_extraction',
            description='Number of successful extractions for each barcode element.',
            plot=bargraph.plot(
                plot_data,
                pconfig={
                    'id': 'barcode_extraction_plot',
                    'title': 'Barcode: Extraction Success',
                    'ylab': 'Number of Reads',
                    'cpswitch_counts_label': 'Number of Reads'
                }
            )
        )
    
    def _add_correction_plot(self):
        """Create barplots showing correction statistics for each barcode type."""
        if not self.correct_data:
            return
            
        # Group data by barcode type
        barcode_groups = defaultdict(dict)
        for s_name, d in self.correct_data.items():
            barcode_name = d.get('barcode_name', 'default')
            counts = d.get('barcode_counts', {})
            correction_details = d.get('correction_details', {})
            
            # Get original sample name
            orig_name = d.get('_sample_name', s_name.rsplit('_', 1)[0])
            
            barcode_groups[barcode_name][orig_name] = {
                'Initial Barcodes': counts.get('unique_barcodes_before', 0),
                'Final Barcodes': counts.get('unique_barcodes_after', 0),
                'Corrected': correction_details.get('barcodes_corrected', 0)
            }
        
        # Create a plot for each barcode type
        for barcode_name, plot_data in barcode_groups.items():
            if not plot_data:  # Skip if no data for this barcode type
                continue
                
            title_prefix = f"{barcode_name} " if barcode_name != 'default' else ''
            
            self.add_section(
                name=f'{title_prefix}Correction Results',
                anchor=f'barcode_correction_{barcode_name}',
                description=f'Summary of barcode correction results for {barcode_name} barcodes.',
                plot=bargraph.plot(
                    plot_data,
                    pconfig={
                        'id': f'barcode_correction_plot_{barcode_name}',
                        'title': f'Barcode: {title_prefix}Correction Results',
                        'ylab': 'Number of Barcodes',
                        'cpswitch_counts_label': 'Number of Barcodes'
                    }
                )
            )
    
    def _add_count_distribution(self):
        """Create line plot showing barcode count distribution."""
        if not self.correct_data:
            return
            
        plot_data = {}
        for s_name, d in self.correct_data.items():
            count_data = d.get('count_data', {})
            hist_data = count_data.get('count_histogram', {})
            if not hist_data:
                continue
                
            # Convert to lists for plotting
            counts = sorted([int(x) for x in hist_data.keys()])
            freqs = [hist_data.get(str(x), 0) for x in counts]
            
            if counts and freqs:  # Only add if we have valid data
                plot_data[s_name] = {
                    'x': counts,
                    'y': freqs
                }
            
        if plot_data:
            self.add_section(
                name='Barcode Frequency Distribution',
                anchor='barcode_distribution',
                description='Distribution of barcode frequencies after correction.',
                plot=linegraph.plot(
                    plot_data,
                    pconfig={
                        'id': 'barcode_freq_plot',
                        'title': 'Barcode: Count Distribution',
                        'xlab': 'Barcode Count',
                        'ylab': 'Number of Barcodes',
                        'xLog': True,
                        'yLog': True
                    }
                )
            ) 