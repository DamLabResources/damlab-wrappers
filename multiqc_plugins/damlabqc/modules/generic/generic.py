from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.generic')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Generic',
            anchor='generic',
            href='',
            info='displays generic metrics and status information.'
        )

        # Initialize data storage
        self.generic_data = dict()
        
        # Find and parse input files
        for f in self.find_log_files('generic'):
            self._parse_generic_log(f)

        # Filter out ignored samples
        self.generic_data = self.ignore_samples(self.generic_data)

        # No samples found
        if len(self.generic_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self.generic_data)} reports")
        
        # Write parsed data to file
        self.write_data_file(self.generic_data, 'multiqc_generic')
        
        # Add to General Statistics table
        self._add_general_stats()
    
    def _parse_generic_log(self, f):
        """Parse generic YAML log file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            
            self.generic_data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            return
    
    def _add_general_stats(self):
        """Add all fields to the general stats table."""
        # Collect all unique keys across all samples
        all_keys = set()
        for data in self.generic_data.values():
            all_keys.update(data.keys())
        
        # Remove sample_name from keys as it's redundant
        if 'sample_name' in all_keys:
            all_keys.remove('sample_name')
        
        # Create headers for each field
        headers = {}
        for key in sorted(all_keys):
            headers[key] = {
                'title': key.replace('_', ' ').title(),
                'description': f'Generic field: {key}',
                'scale': 'RdYlBu',
                'format': '{:}',
                'placement': 1000  # All fields at same priority
            }
            
            # Special formatting for known field types
            if any(status_word in key.lower() for status_word in ['status', 'state']):
                headers[key]['scale'] = 'RdYlGn'
            elif any(numeric_word in key.lower() for numeric_word in ['count', 'total', 'number', 'rate', 'percentage']):
                headers[key]['format'] = '{:,.2f}'
                
        # Add to general stats table
        self.general_stats_addcols(self.generic_data, headers) 