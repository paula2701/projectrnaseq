/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD           } from '../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN           } from '../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_COORD } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FIXMATE       } from '../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_MARKDUP       } from '../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_projectrnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROJECTRNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta          // path: genome fasta
    gtf            // path: genome annotation gtf
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    
    //
    // MODULE: Concatenate FastQ files rom the same sample if required
    //
    ch_samplesheet
        .branch { meta, fastqs ->
            single  : fastqs.size() == 1 || (!meta.single_end && fastqs.size() ==2)
            multiple: true
        }
        .set { ch_fastq }

    CAT_FASTQ (ch_fastq.multiple)
    def ch_reads = CAT_FASTQ.out.reads.mix(ch_fastq.single)

    //
    // MODULE: run fastp (adapter + quality trimming)
    //
    FASTP(
        ch_reads.map { meta, reads -> [ meta, reads, [] ] }, // [] = no adapter fasta -> fastp auto-detects
        false, // discard_trimmed_pass: keep the trimmed reads (we need them for alignment)
        false, // save_trimmed_fail: don't save reads that failed filtering
        false // save_merged: don't merge overlapping PE reads (not wanted for RNA-seq)
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map {_meta, json -> json })

    def ch_trimmed_reads = FASTP.out.reads // input for step 3: alignment

    //
    //MODULE: Prepare HISAT2 reference (splice sites + index), built once
    //
    def ch_gtf = channel.value ([ [id:'genome'], file(gtf, checkIfExists:true) ])
    HISAT2_EXTRACTSPLICESITES(ch_gtf)
    def ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt.first()

    HISAT2_BUILD(
        ch_splicesites.map { meta, ss -> [ meta, file(fasta, checkIfExists:true), file(gtf, checkIfExists:true), ss ] },
        params.hisat2_build_memory
    )
    def ch_hisat2_index = HISAT2_BUILD.out.index.first()

    //
    // MODULE: Align trimmed reads with HISAT2
    //
    HISAT2_ALIGN(ch_trimmed_reads, ch_hisat2_index, ch_splicesites, false)

    //
    // MODULE: Mark duplicates
    //  name sort -> fixmate -m -> coordinate sort -> markdup -> index
    //
    SAMTOOLS_SORT_NAME(HISAT2_ALIGN.out.bam, [[:], [], []], '')         // '' = no index (name-sorted BAM can't be indexed)
    SAMTOOLS_FIXMATE(SAMTOOLS_SORT_NAME.out.bam, [[:], [], []])         // adds ms/MC tags needed by markdup
    SAMTOOLS_SORT_COORD(SAMTOOLS_FIXMATE.out.bam, [[:], [], []], '')    // index only the final BAM
    SAMTOOLS_MARKDUP(SAMTOOLS_SORT_COORD.out.bam, [[:], [], []])
    SAMTOOLS_INDEX(SAMTOOLS_MARKDUP.out.bam)

    def ch_bam_bai = SAMTOOLS_MARKDUP.out.bam.join(SAMTOOLS_INDEX.out.index) // input for next step: quantification

    //
    // MODULE: Mapping stats (now includes duplicate counts)
    //
    SAMTOOLS_FLAGSTAT(ch_bam_bai)

    ch_multiqc_files = ch_multiqc_files
        .mix(HISAT2_ALIGN.out.summary.map { _meta, f -> f })
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { _meta, f -> f })

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file})

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'projectrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'projectrnaseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
