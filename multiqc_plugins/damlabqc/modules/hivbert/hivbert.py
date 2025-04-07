from collections import OrderedDict, defaultdict
import logging
import yaml
import numpy as np

from multiqc import config
from multiqc.plots import heatmap
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.hivbert')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='HIV-BERT',
            anchor='hivbert',
            href='https://huggingface.co/damlab/HIV_BERT',
            info='analyzes HIV sequence predictions using HIV-BERT models.'
        )

        # Initialize data storage
        self.hivbert_data = defaultdict(lambda: {
            'total_sequences': 0,
            'processed_sequences': 0,
            'mean_values': {},
            'high_confidence_counts': {},
            'model_type': None,
            'model_name': None
        })
        
        # Find and parse input files
        for f in self.find_log_files('hivbert'):
            self._parse_hivbert_yaml(f)
            
        # Convert defaultdict to regular dict
        self.hivbert_data = dict(self.hivbert_data)
        
        # Filter out ignored samples
        self.hivbert_data = self.ignore_samples(self.hivbert_data)

        # No samples found
        if len(self.hivbert_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self.hivbert_data)} samples")
        
        # Write parsed data to file
        self.write_data_file(self.hivbert_data, 'multiqc_hivbert')
        
        # Add sections
        self._add_general_stats()
        self._add_prediction_heatmap()
    
    def _parse_hivbert_yaml(self, f):
        """Parse HIV-BERT YAML metrics file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            # Use sample_name from YAML if present, otherwise use filename
            s_name = data.get('sample_name', f['s_name'].split('.')[0])
            s_name = self.clean_s_name(s_name, f)
            
            # Store data
            self.hivbert_data[s_name].update({
                'total_sequences': data['total_sequences'],
                'processed_sequences': data['processed_sequences'],
                'mean_values': data['mean_values'],
                'high_confidence_counts': data['high_confidence_counts'],
                'model_type': data['model_type'],
                'model_name': data['model_name']
            })
            
            log.debug(f"Parsed data for sample {s_name}")
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add sequence counts and prediction rates to general stats table."""
        data = {}
        
        for s_name, d in self.hivbert_data.items():
            data[s_name] = {
                'hivbert-total_sequences': d['total_sequences'],
                'hivbert-processed_sequences': d['processed_sequences'],
                'hivbert-processed_rate': d['processed_sequences'] / d['total_sequences'] if d['total_sequences'] > 0 else 0
            }
            
            # Add mean values for classification models
            if d['model_type'] == 'classification':
                for label, value in d['mean_values'].items():
                    data[s_name][f'hivbert-mean_{label}'] = value
                
                # Add high confidence counts as percentages
                for label, count in d['high_confidence_counts'].items():
                    data[s_name][f'hivbert-high_conf_{label}'] = count / d['processed_sequences'] if d['processed_sequences'] > 0 else 0
        
        # Create headers
        headers = OrderedDict()
        
        headers['hivbert-total_sequences'] = {
            'title': 'HIV-BERT Total',
            'description': 'Total number of input sequences',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        
        headers['hivbert-processed_sequences'] = {
            'title': 'HIV-BERT Processed',
            'description': 'Number of sequences processed by model',
            'format': '{:,.0f}',
            'scale': 'Greens'
        }
        
        headers['hivbert-processed_rate'] = {
            'title': 'HIV-BERT % Processed',
            'description': 'Percentage of sequences successfully processed',
            'max': 1,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.1%}'
        }
        
        # Add headers for classification model stats
        # Get first sample with classification data to determine columns
        class_sample = next((d for d in self.hivbert_data.values() if d['model_type'] == 'classification'), None)
        if class_sample:
            # Add mean prediction headers
            for label in class_sample['mean_values'].keys():
                headers[f'hivbert-mean_{label}'] = {
                    'title': f'Mean {label}',
                    'description': f'Mean prediction probability for {label}',
                    'max': 1,
                    'min': 0,
                    'scale': 'RdYlBu',
                    'format': '{:.2f}'
                }
            
            # Add high confidence headers
            for label in class_sample['high_confidence_counts'].keys():
                headers[f'hivbert-high_conf_{label}'] = {
                    'title': f'High Conf {label}',
                    'description': f'Proportion of sequences with high confidence {label} prediction (≥0.5)',
                    'max': 1,
                    'min': 0,
                    'scale': 'YlOrRd',
                    'format': '{:.1%}'
                }
        
        self.general_stats_addcols(data, headers)
    
    def _add_prediction_heatmap(self):
        """Create heatmap showing predictions across samples."""
        # Collect all prediction labels
        all_labels = set()
        for d in self.hivbert_data.values():
            if d['model_type'] == 'classification':
                all_labels.update(d['mean_values'].keys())
        
        if not all_labels:
            return  # No classification data found
            
        # Prepare heatmap data - ensure proper data structure
        mean_data = {}
        conf_data = {}
        
        for s_name, d in self.hivbert_data.items():
            if d['model_type'] == 'classification':
                # Create row for mean predictions
                mean_data[s_name] = {}
                for label in sorted(all_labels):
                    mean_data[s_name][label] = float(d['mean_values'].get(label, 0))
                
                # Create row for high confidence counts
                conf_data[s_name] = {}
                total = float(d['processed_sequences'])  # Use processed_sequences instead of total
                if total > 0:
                    for label in sorted(all_labels):
                        conf_data[s_name][label] = float(d['high_confidence_counts'].get(label, 0)) / total
        
        if mean_data:
            # Create mean predictions heatmap
            self.add_section(
                name='Mean Predictions',
                anchor='hivbert_heatmap',
                description='Mean prediction probabilities across samples. ' +
                          'Darker colors indicate higher probabilities.',
                helptext='This heatmap shows the average prediction probability for each category across all sequences in each sample.',
                plot=heatmap.plot(
                    mean_data,
                    xcats=sorted(all_labels),
                    ycats=sorted(mean_data.keys()),
                    pconfig={
                        'id': 'hivbert_heatmap',
                        'title': 'HIV-BERT: Mean Predictions',
                        'xTitle': 'Prediction Categories',
                        'yTitle': 'Samples',
                        'min': 0,
                        'max': 1,
                        'square': False,
                        'colstops': [
                            [0, '#FFFFFF'],
                            [0.25, '#FFF7EC'],
                            [0.5, '#FC8D59'],
                            [0.75, '#D73027'],
                            [1, '#4A1486']
                        ],
                        'decimalPlaces': 2,
                        'legend': True,
                        'datalabels': True
                    }
                )
            )
            
            # Create high confidence heatmap
            self.add_section(
                name='High Confidence Predictions',
                anchor='hivbert_confidence',
                description='Proportion of sequences with high confidence predictions (≥0.5) ' +
                          'for each category.',
                helptext='This heatmap shows the proportion of sequences in each sample that have a prediction confidence ≥0.5 for each category.',
                plot=heatmap.plot(
                    conf_data,
                    xcats=sorted(all_labels),
                    ycats=sorted(conf_data.keys()),
                    pconfig={
                        'id': 'hivbert_confidence_heatmap',
                        'title': 'HIV-BERT: High Confidence Predictions',
                        'xTitle': 'Prediction Categories',
                        'yTitle': 'Samples',
                        'min': 0,
                        'max': 1,
                        'square': False,
                        'colstops': [
                            [0, '#FFFFFF'],
                            [0.25, '#FFF7EC'],
                            [0.5, '#FC8D59'],
                            [0.75, '#D73027'],
                            [1, '#4A1486']
                        ],
                        'decimalPlaces': 2,
                        'legend': True,
                        'datalabels': True
                    }
                )
            ) 