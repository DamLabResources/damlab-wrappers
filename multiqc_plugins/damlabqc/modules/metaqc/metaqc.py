from collections import OrderedDict, defaultdict
import logging
import os
import pandas as pd
from typing import Dict, List, Optional, Tuple

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.utils import report

log = logging.getLogger('damlabqc.metaqc')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Meta MultiQC',
            anchor='metaqc',
            href='',
            info='aggregates MultiQC reports across runs by grouping samples.'
        )

        # Initialize data storage
        self.metaqc_data = dict()
        self.all_headers = set()
        
        # Find and load MultiQC general stats files
        for f in self.find_log_files('metaqc'):
            self._parse_stats_file(f)

        # Filter out ignored samples
        self.metaqc_data = self.ignore_samples(self.metaqc_data)

        if len(self.metaqc_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.metaqc_data)} samples")

        # Group samples by directory
        self.grouped_data = self._group_by_directory()
        
        # Write parsed data to file
        self.write_data_file(self.grouped_data, 'multiqc_metaqc')
        
        # Add to General Statistics table
        self._add_general_stats()

    def _parse_stats_file(self, f):
        """Parse MultiQC general stats file."""
        try:
            # Read TSV file
            df = pd.read_csv(f['f'], sep='\t')
            
            # Get directory name as group
            group = os.path.basename(os.path.dirname(os.path.dirname(f['root'])))
            
            # Store data
            for _, row in df.iterrows():
                sample = row['Sample']
                s_name = f"{group}/{sample}"
                s_name = self.clean_s_name(s_name, f)
                
                # Store all data except Sample column
                self.metaqc_data[s_name] = {
                    col: row[col] 
                    for col in df.columns 
                    if col != 'Sample'
                }
                self.all_headers.update(self.metaqc_data[s_name].keys())
                
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            return

    def _group_by_directory(self) -> Dict:
        """Group samples by their directory prefix."""
        groups = defaultdict(dict)
        
        for s_name, data in self.metaqc_data.items():
            group = s_name.split('/')[0]
            sample = s_name.split('/')[1]
            
            # Initialize group if needed
            if group not in groups:
                groups[group] = {
                    header: [] for header in self.all_headers
                }
            
            # Add sample data to group
            for header in self.all_headers:
                if header in data:
                    groups[group][header].append(data[header])
                
        # Calculate aggregates
        aggregated = {}
        for group, data in groups.items():
            aggregated[group] = {}
            for header, values in data.items():
                if not values:
                    continue
                    
                # Calculate mean for numeric values
                try:
                    aggregated[group][header] = sum(values) / len(values)
                except:
                    # For non-numeric, take first value
                    aggregated[group][header] = values[0]
                    
        return aggregated

    def _add_general_stats(self):
        """Add aggregated stats to the general stats table."""
        headers = {}
        for header in sorted(self.all_headers):
            headers[header] = {
                'title': header.replace('_', ' ').title(),
                'description': f'Average {header} across runs',
                'scale': 'RdYlBu',
                'format': '{:,.2f}',
                'shared_key': header.split('-')[0] if '-' in header else None,
            }

        self.general_stats_addcols(
            self.grouped_data,
            headers
        ) 