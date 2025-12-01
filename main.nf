nextflow.enable.dsl=2

// ---------------------------
// Parameters (schema-backed)
// ---------------------------
params.input        = params.input        ?: "${baseDir}/assets/samplesheet.csv"
params.outdir       = params.outdir       ?: "${baseDir}/results"
params.reference    = params.reference    ?: null
params.octopus_args = params.octopus_args ?: ""
params.help         = params.help         ?: false

// ---------------------------
// Help message
// ---------------------------
def helpMessage = """
octopus_call – Small Nextflow pipeline running Octopus on BAM files

Required:
  --input        Path to samplesheet CSV (columns: sample_id,bam)
  --reference    Reference genome FASTA used for Octopus

Optional:
  --outdir       Output directory (default: results)
  --octopus_args Extra arguments passed directly to Octopus

Example:
  nextflow run main.nf \\
    --input assets/samplesheet.csv \\
    --reference data/genome.fa \\
    --outdir results
"""

if (params.help) {
    log.info helpMessage
    System.exit(0)
}

// ---------------------------
// Basic sanity checks
// ---------------------------
if (!params.reference) {
    error "Parameter --reference is required (reference FASTA for Octopus)."
}

if (!file(params.reference).exists()) {
    error "Reference FASTA not found: ${params.reference}"
}

if (!file(params.input).exists()) {
    error "Samplesheet not found: ${params.input}"
}

// ---------------------------
// Channels
// ---------------------------
Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row ->
        if (!row.sample_id) {
            error "Samplesheet is missing 'sample_id' column or has empty values."
        }
        if (!row.bam) {
            error "Samplesheet is missing 'bam' column or has empty values."
        }
        def sample_id = row.sample_id as String
        def bam_path  = file(row.bam)

        if (!bam_path.exists()) {
            error "BAM file for sample '${sample_id}' not found: ${bam_path}"
        }

        tuple(sample_id, bam_path)
    }
    .set { ch_samples }

// ---------------------------
// Processes
// ---------------------------
process OCTOPUS {

    tag "${sample_id}"

    publishDir "${params.outdir}/octopus", mode: 'copy'

    cpus 4
    memory '8 GB'

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}.octopus.vcf.gz", emit: vcf
    path "${sample_id}.octopus.vcf.gz.tbi", optional: true, emit: tbi

    when:
    bam

    script:
    """
    octopus \\
      -R ${params.reference} \\
      -I ${bam} \\
      -o ${sample_id}.octopus.vcf.gz \\
      ${params.octopus_args}
    """
}

// ---------------------------
// Workflow definition
// ---------------------------
workflow {

    log.info "Launching octopus_call pipeline"
    log.info "Samplesheet : ${params.input}"
    log.info "Reference   : ${params.reference}"
    log.info "Outdir      : ${params.outdir}"

    OCTOPUS(ch_samples)
}
