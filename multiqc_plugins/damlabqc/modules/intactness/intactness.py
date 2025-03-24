from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.plots import bargraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.intactness')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Intactness',
            anchor='intactness',
            href='',
            info='analyzes HIV sequence intactness metrics.'
        )

        # Initialize data storage
        self.intactness_data = dict()
        
        # Find and parse logs
        for f in self.find_log_files('intactness'):
            self._parse_intactness_log(f)

        # Filter out ignored samples
        self.intactness_data = self.ignore_samples(self.intactness_data)

        # No samples found
        if len(self.intactness_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self.intactness_data)} reports")
        
        # Write parsed data to file
        self.write_data_file(self.intactness_data, 'multiqc_intactness')

        # Add to general stats table
        self._add_general_stats()
        
        # Add barplot section
        self._add_barplot()
    
    def _parse_intactness_log(self, f):
        """Parse intactness YAML file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            
            self.intactness_data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add intactness metrics to the general stats table."""
        general_stats_data = {}
        
        for s_name, data in self.intactness_data.items():
            general_stats_data[s_name] = {
                'total_reads': data['total_reads'],
                'countable_reads': data['countable_reads'],
                'intact_reads': data['intact_reads'],
                'percent_countable': data['percent_countable'],
                'percent_intact': data['percent_intact']
            }

        headers = OrderedDict()
        headers['percent_intact'] = {
            'title': '% Intact',
            'description': 'Percentage of countable reads that are intact',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:,.1f}%',
            'placement': 2400
        }
        headers['total_reads'] = {
            'title': 'Total Reads',
            'description': 'Total number of reads',
            'scale': 'GnBu',
            'format': '{:,.0f}',
            'placement': 2000,
            'hidden': True
        }
        headers['countable_reads'] = {
            'title': 'Countable',
            'description': 'Number of reads within countable length range',
            'scale': 'YlGn',
            'format': '{:,.0f}',
            'placement': 2100,
            'hidden': True
        }
        headers['intact_reads'] = {
            'title': 'Intact',
            'description': 'Number of reads within intact length range',
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'placement': 2200,
            'hidden': True
        }
        headers['percent_countable'] = {
            'title': '% Countable',
            'description': 'Percentage of total reads that are countable',
            'max': 100,
            'min': 0,
            'scale': 'YlOrRd',
            'format': '{:,.1f}%',
            'placement': 2300,
            'hidden': True
        }
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_barplot(self):
        """Add stacked barplot showing intact vs non-intact reads."""
        plot_data = {}
        cats = ['Intact', 'Not Intact', 'Not Countable']
        
        for s_name, d in self.intactness_data.items():
            plot_data[s_name] = dict()
            plot_data[s_name]['Intact'] = d['intact_reads']
            plot_data[s_name]['Not Intact'] = d['countable_reads'] - d['intact_reads']
            plot_data[s_name]['Not Countable'] = d['total_reads'] - d['countable_reads']

        self.add_section(
            name='Read Intactness',
            anchor='intactness-distribution',
            description='Distribution of reads by intactness classification.',
            helptext='''
            This plot shows the distribution of reads across three categories:
            * **Intact**: Reads within both countable and intact length ranges
            * **Not Intact**: Reads within countable but not intact length range
            * **Not Countable**: Reads outside the countable length range
            
            Use the buttons above the plot to switch between showing absolute read counts and percentages.
            ''',
            plot=bargraph.plot(
                plot_data,
                cats,
                pconfig={
                    'id': 'intactness_distribution',
                    'title': 'Intactness: Read Distribution',
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