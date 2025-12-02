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
  --input        Path to samplesheet CSV (columns: sample_id,bam[,bai])
  --reference    Reference genome FASTA (.fa) used for Octopus

Optional:
  --outdir       Output directory (default: results)
  --octopus_args Extra arguments passed directly to Octopus

Example:
  nextflow run main.nf \\
    --input assets/samplesheet.csv \\
    --reference data/GRCh38_chr22.fa \\
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

if (!file(params.reference + ".fai").exists()) {
    error "Reference FASTA index (.fai) not found: ${params.reference}.fai"
}

if (!file(params.input).exists()) {
    error "Samplesheet not found: ${params.input}"
}

// Create file objects
def reference_fa  = file(params.reference)
def reference_fai = file(params.reference + ".fai")

// Create channel for reference + fai
Channel.value(reference_fa).set { ch_reference_fa }
Channel.value(reference_fai).set { ch_reference_fai }

// ---------------------------
// Channels
// ---------------------------
Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row ->
        if (!row.sample_id) error "Samplesheet missing 'sample_id'"
        if (!row.bam)       error "Samplesheet missing 'bam'"

        def sample_id = row.sample_id as String
        def bam       = file(row.bam)

        if (!bam.exists()) {
            error "BAM not found for ${sample_id}: ${bam}"
        }

        // Use bai column if present, else try <bam>.bai
        def bai = row.bai ? file(row.bai) : file("${bam}.bai")

        if (!bai.exists()) {
            // fallback: <bam_base>.bai
            def altBai = file("${bam.parent}/${bam.baseName}.bai")
            if (altBai.exists()) bai = altBai
            else error "BAI not found for ${sample_id}. Tried: ${bai} and ${altBai}"
        }

        tuple(sample_id, bam, bai)
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
    tuple val(sample_id), path(bam), path(bai)
    path reference_fa
    path reference_fai

    output:
    path "${sample_id}.octopus.vcf.gz", emit: vcf
    path "${sample_id}.octopus.vcf.gz.tbi", optional: true, emit: tbi

    script:
    """
    echo "=== Running Octopus on ${sample_id} ==="
    echo "BAM       : ${bam}"
    echo "BAI       : ${bai}"
    echo "Reference : ${reference_fa}"
    echo "Index     : ${reference_fai}"

    # Ensure BAM and index have matching names inside work dir
    if [ ! -e "${bam}.bai" ]; then
      ln -s "${bai}" "${bam}.bai"
    fi

    octopus \\
      -R ${reference_fa} \\
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

    OCTOPUS(ch_samples, ch_reference_fa, ch_reference_fai)
}
