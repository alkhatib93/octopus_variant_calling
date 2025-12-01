# octopus_call

Small Nextflow pipeline that runs the [Octopus](https://github.com/luntergroup/octopus) variant caller on one or more BAM files, using a simple samplesheet as input.

- **Input**: CSV samplesheet with at least `sample_id` and `bam` columns
- **Output**: One gzipped VCF (`*.octopus.vcf.gz`) per sample + optional index

## Features

- Minimal, focused wrapper around Octopus
- Samplesheet-based input (easy to integrate with other tooling)
- Schema-backed parameters via `nextflow_schema.json`
- Simple, extensible structure (add annotation or filtering downstream)

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Octopus in `$PATH`, or use a container/Conda env via `nextflow.config`
- POSIX-like environment (Linux recommended)

## Quickstart

1. Prepare a reference FASTA (and indexes) for Octopus:

   ```bash
   # example
   data/genome.fa
