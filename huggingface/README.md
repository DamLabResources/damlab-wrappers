# HuggingFace Wrappers

This directory contains wrappers for HuggingFace models.

## HIV-BERT Wrapper

The HIV-BERT wrapper allows you to run sequences through HuggingFace HIV-BERT models.

### Environment Setup

The wrapper supports two ways to set up the environment:

1. **Default Environment**: Uses the latest versions of PyTorch and transformers.
2. **Custom Environment**: Allows you to specify a custom environment with pinned versions.

#### Using the Default Environment

By default, the wrapper will use the environment specified in `hiv-bert/environment.yaml`.

#### Using a Custom Environment

For clusters with specific CUDA requirements (like Picotte), you can create a custom environment:

1. Build the Picotte-specific environment:
   ```bash
   cd /path/to/damlab-wrappers/huggingface
   make install-picotte
   ```

2. Configure your Snakemake workflow to use the custom environment:
   ```yaml
   # In your config.yaml
   hiv_bert_env: /path/to/damlab-wrappers/huggingface/venv/picotte
   ```

3. Run your workflow:
   ```bash
   snakemake --use-conda
   ```

### Available Models

- `damlab/hiv_bert`: Base HIV-BERT model for embeddings
- `damlab/HIV_BERT`: Alternative name for the base model
- `damlab/HIV_V3_bodysite`: Predicts HIV V3 sequence body site origin
- `damlab/HIV_V3_coreceptor`: Predicts HIV V3 sequence co-receptor tropism

### Parameters

- `model_name`: Name of the HuggingFace model to use
- `model_directory`: Directory to cache models
- `min_length`: Minimum sequence length
- `max_length`: Maximum sequence length
- `force_cpu`: Force CPU mode (set to True to bypass GPU)
- `custom_env`: Path to custom environment

### Example Usage

```python
rule HIV_V3_bodysite:
    input:
        'sliced/{sample}.V3.fastq',
    output:
        'V3/{sample}.bodysite.csv',
        metrics =  'V3/{sample}.bodysite.yaml',
    params:
        model_name = 'damlab/HIV_V3_bodysite',
        model_directory = 'model_cache',
        min_length = 30,
        max_length = 36,
        sample_name = lambda wildcards: wildcards['sample'],
        force_cpu = False,
        custom_env = config.get('hiv_bert_env', None)
    threads: max(int(workflow.cores/2),1)
    wrapper:
        f"file:{DL_PREFIX}/huggingface/hiv-bert"
``` 