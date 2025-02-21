from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.plots import bargraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.deletion_frequency')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Deletion Frequency',
            anchor='deletion_frequency',
            href='',
            info='analyzes deletion frequencies from BAM files.'
        )

        # Initialize data storage
        self._deletion_freq_data = dict()
        self._regions = set()
        
        self._collect_log_files()
        self._add_general_stats()
        self._add_sections()
    
    def _collect_log_files(self):
        """Collect data from log files found by multiqc."""
        for f in self.find_log_files('deletion_frequency'):
            self._parse_deletion_freq_log(f)

        # Filter out ignored samples
        self._deletion_freq_data = self.ignore_samples(self._deletion_freq_data)

        # No samples found
        if len(self._deletion_freq_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self._deletion_freq_data)} reports")
        
        # Write parsed report data
        self.write_data_file(self._deletion_freq_data, 'multiqc_deletion_frequency')
    
    def _parse_deletion_freq_log(self, f):
        """Parse single deletion frequency log file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            region = data.get('region_name', 'region')
            
            # Create nested structure if sample doesn't exist
            if s_name not in self._deletion_freq_data:
                self._deletion_freq_data[s_name] = {}
            
            # Store data under sample and region
            self._deletion_freq_data[s_name][region] = data
            self._regions.add(region)
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add to the general stats table at the top of the report."""
        general_stats_data = {}
        headers = OrderedDict()
        
        # Create headers for each region
        for region in sorted(self._regions):
            headers[f'deletion_frequency_{region}'] = {
                'title': f'{region} Del. Freq.',
                'description': f'Deletion frequency for {region}',
                'max': 1,
                'min': 0,
                'scale': 'RdYlBu-rev',
                'format': '{:,.2%}',
                'placement': 1000
            }
            headers[f'reads_with_deletion_{region}'] = {
                'title': f'{region} Reads w/Del',
                'description': f'Number of reads containing the deletion in {region}',
                'scale': 'PuRd',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1100
            }
            headers[f'reads_covering_required_{region}'] = {
                'title': f'{region} Reads Covering',
                'description': f'Number of reads covering the required region in {region}',
                'scale': 'YlGn',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1200
            }
            headers[f'total_reads_{region}'] = {
                'title': f'{region} Total Reads',
                'description': f'Total number of reads in {region}',
                'scale': 'GnBu',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1300
            }

        # Fill in the data
        for s_name, regions in self._deletion_freq_data.items():
            general_stats_data[s_name] = {}
            for region, data in regions.items():
                general_stats_data[s_name].update({
                    f'total_reads_{region}': data['total_reads'],
                    f'reads_covering_required_{region}': data['reads_covering_required'],
                    f'reads_with_deletion_{region}': data['reads_with_deletion'],
                    f'deletion_frequency_{region}': data['deletion_frequency']
                })
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_sections(self):
        """Add sections to MultiQC report."""
        for region in sorted(self._regions):
            self._add_deletion_stats(region)
    
    def _add_deletion_stats(self, region):
        """Add deletion statistics plot for a specific region."""
        deletion_data = {}
        
        for s_name, regions in self._deletion_freq_data.items():
            if region in regions:
                d = regions[region]
                deletion_data[s_name] = {
                    'With Deletion': d['reads_with_deletion'],
                    'Without Deletion': d['reads_covering_required'] - d['reads_with_deletion'],
                    'Not Covering': d['total_reads'] - d['reads_covering_required']
                }

        self.add_section(
            name=f'Deletion Statistics - {region}',
            anchor=f'deletion-frequency-stats-{region}',
            description=f'Distribution of reads with respect to deletion and coverage for {region}',
            plot=bargraph.plot(deletion_data,
                             pconfig={
                                 'id': f'deletion_frequency_stats_{region}',
                                 'title': f'Deletion Frequency: Read Statistics - {region}',
                                 'ylab': 'Number of Reads',
                                 'stacking': 'normal'
                             })
        ) 