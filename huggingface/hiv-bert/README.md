# HIV-BERT Wrapper

A Snakemake wrapper for running sequences through HIV-BERT models from HuggingFace. This wrapper handles both DNA and amino acid sequences, automatically detecting and processing them appropriately.

## Models

The wrapper supports the following HIV-BERT models:

- `damlab/hiv_bert`: Base HIV-BERT model that provides sequence embeddings
- `damlab/HIV_V3_bodysite`: Predicts HIV V3 sequence body site origin
- `damlab/HIV_V3_coreceptor`: Predicts HIV V3 sequence co-receptor tropism

## Input

The wrapper accepts the following input formats:
- FASTA files (`.fa`, `.fasta`)
- FASTQ files (`.fq`, `.fastq`)
- BAM/SAM files (`.bam`, `.sam`)

Both DNA and amino acid sequences are supported. DNA sequences will be automatically translated to amino acids before processing.

## Output

The output is a CSV file containing:
- For embedding models: sequence IDs and their 768-dimensional embeddings
- For classification models: sequence IDs and probabilities for each class

## Parameters

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| model_name | HuggingFace model name | string | 'damlab/hiv_bert' |
| model_directory | Directory to cache models | string | None |
| mapped_only | Only process mapped reads (BAM/SAM) | boolean | False |
| translate_frame | Reading frame for DNA translation | int | 0 |
| min_length | Minimum amino acid sequence length | int | 10 |
| max_length | Maximum amino acid sequence length | int | 256 |
| batch_size | Batch size for processing sequences | int | 32 |

## Example Usage

### Basic Usage

```python
rule predict_bodysite:
    input:
        "sequences.fasta"  # Can be DNA or AA sequences
    output:
        "predictions.csv"
    params:
        model_name = "damlab/HIV_V3_bodysite"
    wrapper:
        "file:path/to/damlab-wrappers/huggingface/hiv-bert"
```

### With Model Caching

```python
rule get_embeddings:
    input:
        "sequences.fasta"
    output:
        "embeddings.csv"
    params:
        model_name = "damlab/hiv_bert",
        model_directory = "model_cache",  # Cache models here
        min_length = 20,
        max_length = 512
    wrapper:
        "file:path/to/damlab-wrappers/huggingface/hiv-bert"
```

### Processing BAM Files

```python
rule process_alignments:
    input:
        "aligned.bam"
    output:
        "predictions.csv"
    params:
        model_name = "damlab/HIV_V3_coreceptor",
        mapped_only = True,
        translate_frame = 0
    wrapper:
        "file:path/to/damlab-wrappers/huggingface/hiv-bert"
```

## Output Format

### Embedding Model Output
```csv
id,dim_0,dim_1,...,dim_767
seq1,0.123,0.456,...,0.789
seq2,-0.234,0.567,...,0.890
```

### Bodysite Model Output
```csv
id,periphery-tcell,periphery-monocyte,CNS,breast-milk,female-genitals,male-genitals,gastric,lung,organ
seq1,0.92,0.03,0.01,0.00,0.02,0.01,0.00,0.00,0.01
seq2,0.01,0.85,0.02,0.00,0.05,0.02,0.03,0.01,0.01
```

### Coreceptor Model Output
```csv
id,CCR5,CXCR4,CCR5+CXCR4,neither
seq1,0.95,0.02,0.02,0.01
seq2,0.03,0.90,0.05,0.02
```

## Requirements

The wrapper requires the following dependencies:
- Python ≥3.10
- pandas
- pysam
- biopython
- transformers
- torch
- sentencepiece
- protobuf

These dependencies are automatically managed through the conda environment specified in the wrapper.

## Notes

1. DNA sequences are automatically detected and translated to amino acids
2. Case is ignored for both DNA and amino acid sequences
3. Invalid characters in amino acid sequences are replaced with 'X'
4. Short sequences (below min_length) are filtered out
5. Long sequences (above max_length) are truncated
6. Models are cached if model_directory is specified
7. GPU is used automatically if available

## Author

Will Dampier (wnd22@drexel.edu)

## License

MIT

## Building the Environment on GPU Nodes

This wrapper requires PyTorch with CUDA support, which can only be built on nodes with GPU access. The included makefiles use `srun` to build the environment on a GPU node.

### Main Wrapper

To build the main wrapper environment:

```bash
cd damlab-wrappers/huggingface/hiv-bert
make
```

This will create a virtual environment in the `venv` directory using the `environment.picotte.yaml` file.

### Test Environment

To build and run the tests:

```bash
cd damlab-wrappers/huggingface/hiv-bert/test
make
make test
```

This will create a virtual environment in the `venv` directory and run the tests using the GPU.

### Cleaning Up

To clean up the environments:

```bash
# For the main wrapper
cd damlab-wrappers/huggingface/hiv-bert
make clean

# For the tests
cd damlab-wrappers/huggingface/hiv-bert/test
make clean
``` 