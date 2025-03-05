from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.plots import bargraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.primercheck')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Primercheck',
            anchor='primercheck',
            href='',
            info='analyzes primer amplification metrics.'
        )

        # Initialize data storage
        self.primercheck_data = dict()
        
        # Find and parse input files
        for f in self.find_log_files('primercheck'):
            self._parse_primercheck_yaml(f)
            
        # Filter out ignored samples
        self.primercheck_data = self.ignore_samples(self.primercheck_data)

        # No samples found
        if len(self.primercheck_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self.primercheck_data)} reports")
        
        # Write parsed data to file
        self.write_data_file(self.primercheck_data, 'multiqc_primercheck')
        
        # Add sections
        self._add_general_stats()
        self._add_pairwise_plot()
    
    def _parse_primercheck_yaml(self, f):
        """Parse primercheck YAML file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            # Use sample_name from YAML if present, otherwise use filename
            s_name = data.get('sample_name', f['s_name'].split('.')[0])
            s_name = self.clean_s_name(s_name, f)
            
            self.primercheck_data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add primer hit frequencies to general stats table."""
        data = {}
        
        for s_name, d in self.primercheck_data.items():
            # Calculate hit rates for each primer
            total_seqs = d['total_sequences']
            if total_seqs > 0:  # Avoid division by zero
                data[s_name] = {
                    f"primer_{name}": hits/total_seqs 
                    for name, hits in d['primer_hits'].items()
                }
                # Add total sequences
                data[s_name]['total_sequences'] = total_seqs
        
        # Create headers
        headers = OrderedDict()
        
        # Add total sequences first
        headers['total_sequences'] = {
            'title': 'Total Seqs',
            'description': 'Total number of sequences analyzed',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        
        # Add a column for each primer
        for primer in self.primercheck_data[list(self.primercheck_data.keys())[0]]['primer_hits'].keys():
            headers[f'primer_{primer}'] = {
                'title': f'{primer}',
                'description': f'Hit rate for primer {primer}',
                'max': 1,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:.1%}'
            }
            
        self.general_stats_addcols(data, headers)
    
    def _add_pairwise_plot(self):
        """Create barplot showing pairwise primer hits."""
        plot_data = {}
        
        for s_name, d in self.primercheck_data.items():
            plot_data[s_name] = d['pairwise_hits']
            
        if plot_data:
            self.add_section(
                name='Pairwise Primer Hits',
                anchor='primercheck_pairwise',
                description='Number of sequences that match each pair of primers.',
                plot=bargraph.plot(
                    plot_data,
                    pconfig={
                        'id': 'primercheck_pairwise_plot',
                        'title': 'Primercheck: Pairwise Primer Hits',
                        'ylab': 'Number of Sequences',
                        'cpswitch_counts_label': 'Number of Sequences'
                    }
                )
            ) 