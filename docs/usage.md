
---

## `docs/USAGE.md`

```markdown
# Usage

This document describes how to run the **octopus_call** pipeline and how to structure the input samplesheet.

---

## Samplesheet format

The pipeline expects a CSV file with at least the following columns:

- `sample_id` – Unique identifier for the sample
- `bam`       – Path to the BAM file for that sample

Example (`assets/samplesheet.csv`):

```csv
sample_id,bam
DEMO_01,/data/bams/DEMO_01.bam
DEMO_02,/data/bams/DEMO_02.bam
