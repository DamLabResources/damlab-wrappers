from collections import OrderedDict, defaultdict
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

        # Initialize data storage with defaultdict for merging
        self.primercheck_data = defaultdict(lambda: {
            'total_sequences': 0,
            'primer_hits': defaultdict(int),
            'pairwise_hits': defaultdict(int)
        })
        
        # Find and parse input files
        for f in self.find_log_files('primercheck'):
            self._parse_primercheck_yaml(f)
            
        # Convert defaultdict to regular dict for downstream processing
        self.primercheck_data = dict(self.primercheck_data)
        
        # Filter out ignored samples
        self.primercheck_data = self.ignore_samples(self.primercheck_data)

        # No samples found
        if len(self.primercheck_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self.primercheck_data)} samples")
        
        # Write parsed data to file
        self.write_data_file(self.primercheck_data, 'multiqc_primercheck')
        
        # Add sections
        self._add_general_stats()
        self._add_pairwise_plot()
    
    def _parse_primercheck_yaml(self, f):
        """Parse primercheck YAML file and merge with existing data if needed."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            # Use sample_name from YAML if present, otherwise use filename
            s_name = data.get('sample_name', f['s_name'].split('.')[0])
            s_name = self.clean_s_name(s_name, f)
            
            # Merge data with existing sample data
            sample_data = self.primercheck_data[s_name]
            
            # Add total sequences
            sample_data['total_sequences'] += data['total_sequences']
            
            # Merge primer hits
            for primer, hits in data['primer_hits'].items():
                sample_data['primer_hits'][primer] += hits
            
            # Merge pairwise hits
            for pair, hits in data['pairwise_hits'].items():
                sample_data['pairwise_hits'][pair] += hits
            
            log.debug(f"Merged data for sample {s_name}")
            
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
        
        # Get all unique primers across all samples
        all_primers = set()
        for d in self.primercheck_data.values():
            all_primers.update(d['primer_hits'].keys())
        
        # Add a column for each primer
        for primer in sorted(all_primers):
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