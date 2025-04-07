from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.plots import bargraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.slice')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Region Slice',
            anchor='slice',
            href='',
            info='analyzes metrics from region-sliced BAM files.'
        )

        # Initialize data storage
        self.slice_data = dict()
        self.regions = set()
        self.region_data = {}
        self.samples = set()
        
        # Find and parse logs
        for f in self.find_log_files('slice'):
            self._parse_slice_log(f)

        # Filter out ignored samples
        self.region_data = self.ignore_samples(self.region_data)

        # No samples found
        if len(self.region_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found reports for {len(self.samples)} samples across {len(self.regions)} regions")
        
        # Write parsed data to file
        self.write_data_file(self.slice_data, 'multiqc_slice')

        # Add to general stats table
        self._add_general_stats()
        
        # Add barplot section
        self._add_barplots()
    
    def _parse_slice_log(self, f):
        """Parse slice YAML metrics file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
            
            # Use sample_name from YAML if available, otherwise use filename
            s_name = data.get('sample_name', f['s_name'])
            if not s_name:
                s_name = f['s_name']
                
            s_name = self.clean_s_name(s_name, f)
            
            # Store region name for later use
            region_name = data.get('region_name', data.get('region', 'unknown'))
            self.regions.add(region_name)
            self.samples.add(s_name)
            
            # Create a unique key for this sample/region combo
            key = f"{s_name}|{region_name}"
            
            # Add region name to the data for later reference
            data['_region_name'] = region_name
            data['_sample_name'] = s_name
            
            # Store data with both methods of access
            self.slice_data[key] = data
            
            # Group by region for stats display
            if region_name not in self.region_data:
                self.region_data[region_name] = {}
            self.region_data[region_name][s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            log.debug(f"Exception details: {str(e)}")
    
    def _add_general_stats(self):
        """Add slice metrics to the general stats table."""
        general_stats_data = {}
        
        # Initialize data structure for each sample
        for s_name in self.samples:
            general_stats_data[s_name] = {}
        
        # For each region, add columns with region name appended
        for region_name, samples in self.region_data.items():
            for s_name, data in samples.items():
                # Make sure sample exists in general_stats_data
                if s_name not in general_stats_data:
                    general_stats_data[s_name] = {}
                    
                # Add data with region name in column name
                general_stats_data[s_name][f'total_segments_{region_name}'] = data['total_segments_processed']
                general_stats_data[s_name][f'overlapping_{region_name}'] = data['segments_overlapping_region']
                
                # Calculate percentage
                if data['total_segments_processed'] > 0:
                    pct = (data['segments_overlapping_region'] / data['total_segments_processed'] * 100)
                else:
                    pct = 0
                    
                general_stats_data[s_name][f'pct_overlapping_{region_name}'] = pct
                general_stats_data[s_name][f'region_{region_name}'] = data.get('region', '')
                general_stats_data[s_name][f'used_index_{region_name}'] = data.get('used_index', False)

        # Create headers for each region
        headers = OrderedDict()
        
        # Sort regions for consistent display
        sorted_regions = sorted(self.regions)
        
        for region_name in sorted_regions:
            headers[f'pct_overlapping_{region_name}'] = {
                'title': f'% Ovlp ({region_name})',
                'description': f'Percentage of processed segments that overlap the {region_name} region',
                'max': 100,
                'min': 0,
                'scale': 'YlGnBu',
                'format': '{:,.1f}%',
                'placement': 1100 + (sorted_regions.index(region_name) * 10)
            }
            headers[f'overlapping_{region_name}'] = {
                'title': f'Reads ({region_name})',
                'description': f'Number of segments overlapping the {region_name} region',
                'scale': 'PuBu',
                'format': '{:,.0f}',
                'placement': 1200 + (sorted_regions.index(region_name) * 10),
                'hidden': len(self.regions) > 2  # Hide if there are many regions
            }
            headers[f'total_segments_{region_name}'] = {
                'title': f'Total ({region_name})',
                'description': f'Total number of segments processed for {region_name} region',
                'scale': 'Blues',
                'format': '{:,.0f}',
                'placement': 1300 + (sorted_regions.index(region_name) * 10),
                'hidden': True
            }
            headers[f'region_{region_name}'] = {
                'title': f'Region ({region_name})',
                'description': f'Target genomic region for {region_name}',
                'placement': 1400 + (sorted_regions.index(region_name) * 10),
                'hidden': True
            }
            headers[f'used_index_{region_name}'] = {
                'title': f'Used Index ({region_name})',
                'description': f'Whether a BAM index was used for fetching {region_name} region',
                'placement': 1500 + (sorted_regions.index(region_name) * 10),
                'hidden': True
            }
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_barplots(self):
        """Add barplots showing read statistics."""
        # 1. Create combined datasets for each region
        for region_name in sorted(self.regions):
            plot_data = {}
            cats = ['Overlapping Region', 'Outside Region']
            
            # Get data for this region across all samples
            for s_name, data in self.region_data.get(region_name, {}).items():
                plot_data[s_name] = dict()
                plot_data[s_name]['Overlapping Region'] = data['segments_overlapping_region']
                plot_data[s_name]['Outside Region'] = data['total_segments_processed'] - data['segments_overlapping_region']

            if plot_data:  # Only create plot if we have data
                self.add_section(
                    name=f'Region: {region_name}',
                    anchor=f'slice-region-{region_name}',
                    description=f'Distribution of reads overlapping the {region_name} region.',
                    helptext=f'''
                    This plot shows the number of reads that overlap the {region_name} region.
                    * **Overlapping Region**: Reads that overlap the specified region
                    * **Outside Region**: Reads that don't overlap the specified region
                    
                    Use the buttons above the plot to switch between showing absolute read counts and percentages.
                    ''',
                    plot=bargraph.plot(
                        plot_data,
                        cats,
                        pconfig={
                            'id': f'slice_region_{region_name}',
                            'title': f'Slice: {region_name} Region',
                            'ylab': 'Number of Reads',
                            'stacking': 'normal',
                            'cpswitch': True,
                            'cpswitch_counts_label': 'Number of Reads',
                            'cpswitch_percent_label': 'Percentage',
                            'tt_decimals': 0,
                            'tt_percentages': True
                        }
                    )
                )
        
        # 2. Create a summary plot if we have multiple regions
        if len(self.regions) > 1:
            # Prepare data for a summary plot with region percentages
            summary_data = {}
            
            for s_name in self.samples:
                summary_data[s_name] = {}
                
                for region_name in sorted(self.regions):
                    if s_name in self.region_data.get(region_name, {}):
                        data = self.region_data[region_name][s_name]
                        if data['total_segments_processed'] > 0:
                            pct = (data['segments_overlapping_region'] / data['total_segments_processed'] * 100)
                        else:
                            pct = 0
                        summary_data[s_name][region_name] = pct
                    else:
                        summary_data[s_name][region_name] = 0
            
            # Add a summary section comparing all regions
            if summary_data:
                self.add_section(
                    name='Region Comparison',
                    anchor='slice-region-comparison',
                    description='Comparison of overlap percentages across all regions.',
                    helptext='''
                    This plot shows the percentage of reads that overlap each targeted genomic region.
                    Higher percentages indicate more reads mapping to that specific region.
                    ''',
                    plot=bargraph.plot(
                        summary_data,
                        sorted(self.regions),
                        pconfig={
                            'id': 'slice_region_comparison',
                            'title': 'Slice: Region Comparison',
                            'ylab': 'Percentage of Reads Overlapping',
                            'cpswitch_c_active': False,
                            'ymax': 100,
                            'tt_suffix': '%',
                            'tt_decimals': 1
                        }
                    )
                ) 