/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// 20240419 s

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// https://github.com/nextflow-io/nf-validation

// Print help message, supply typical command line usage for the pipeline 
if (params.help){
  log.info paramsHelp("nextflow run my_pipeline --input input.csv")
  exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// Create a new channel of metadata from a sample sheet
// NB: `input` corresponds to `params.input` and associated sample sheet schema
ch_input = Channel.fromSamplesheet("input")


// 20240419 e



include { SOMALIER_RELATE } from '../modules/nf-core/somalier/relate/main'
include { VERIFYBAMID_VERIFYBAMID } from '../modules/nf-core/verifybamid/verifybamid/main'
include { VERIFYBAMID_VERIFYBAMID2 } from '../modules/nf-core/verifybamid/verifybamid2/main'
include { CUTADAPT 					} from '../modules/nf-core/cutadapt/main'
include { FASTP as FASTP_TRIM		} from '../modules/nf-core/fastp/main'
include { FASTP as FASTP_RAW 		} from '../modules/nf-core/fastp/main'


include { FASTQ_ALIGN_BWA } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { GATK4SPARK_MARKDUPLICATES } from '../modules/nf-core/gatk4spark/markduplicates/main'
include { GATK4_MARKDUPLICATES } from '../modules/nf-core/gatk4/markduplicates/main'


include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nextbigthing_pipeline'


include { FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON } from '../subworkflows/local/fastq_align_bwamem_mem2_dragmap_sentieon'
include { BAM_MERGE_INDEX_SAMTOOLS                 } from '../subworkflows/local/bam_merge_index_samtools/main'
include { BAM_MARKDUPLICATES                       } from '../subworkflows/local/bam_markduplicates/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING     } from '../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL       } from '../modules/nf-core/samtools/convert/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NEXTBIGTHING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

/*	
	params.step       = 'mapping'
	params.aligner    = "bwa-mem"
	params.fasta      = "/ess/dlstibm/Database/9606.Homo_Sapiens/GRCh38/Custom/sarek.GATK.GRCh38.aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
	params.fasta_fai  = "/ess/dlstibm/Database/9606.Homo_Sapiens/GRCh38/Custom/sarek.GATK.GRCh38.aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai"
	params.dict       = "/ess/dlstibm/Database/9606.Homo_Sapiens/GRCh38/Custom/sarek.GATK.GRCh38.aws/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict" 
	params.bwa        = "/ess/dlstibm/Database/9606.Homo_Sapiens/GRCh38/Custom/sarek.GATK.GRCh38.aws/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex"
*/	
	
	
	// Initialize fasta file with meta map:
	fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

    // Gather built indices or get them from the params
    // Built from the fasta file:
	dict        = params.dict       ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect() : Channel.empty() 
	fasta_fai   = params.fasta_fai  ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it ] }.collect() : Channel.empty()
	bwa         = params.bwa        ? Channel.fromPath(params.bwa).map{ it -> [ [id:'bwa'], it ] }.collect() : Channel.empty()
//	bwamem2     = params.bwamem2    ? Channel.fromPath(params.bwamem2).map{ it -> [ [id:'bwamem2'], it ] }.collect() : PREPARE_GENOME.out.bwamem2
//	dragmap     = params.dragmap    ? Channel.fromPath(params.dragmap).map{ it -> [ [id:'dragmap'], it ] }.collect() : PREPARE_GENOME.out.hashtable


	aligner = params.aligner 

    // Gather index for mapping given the chosen aligner
	index_alignment = (aligner == "bwa-mem" || aligner == "sentieon-bwamem") ? bwa :
		        aligner == "bwa-mem2" ? bwamem2 :
				        dragmap


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    CUTADAPT (
        ch_samplesheet
    )

    // CUTADAPT.out.log.view()
    // CUTADAPT.out.view() : ERR
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]})
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

	ch_samplesheet.view()
		
	
	
	discard_trimmed_pass = false
	save_trimmed_fail = false
	save_merged = false

	FASTP_RAW (
			ch_samplesheet,
			[],
			discard_trimmed_pass,
			save_trimmed_fail,
			save_merged
		)

	    ch_multiqc_files = ch_multiqc_files.mix(FASTP_RAW.out.json.collect{ meta, json -> json })
	    ch_multiqc_files = ch_multiqc_files.mix(FASTP_RAW.out.html.collect{ meta, html -> html })
	    ch_versions = ch_versions.mix(FASTP_RAW.out.versions.first())

	
	if( params.trim_fastq || params.split_fastq > 0 ){

		FASTP_TRIM (
			CUTADAPT.out.reads,
			[],
			discard_trimmed_pass,
			save_trimmed_fail,
			save_merged
		)

		if (params.split_fastq) {

        	reads_for_alignment = FASTP_TRIM.out.reads.map{ meta, reads ->
				read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
				[ meta + [ n_fastq: read_files.size() ], read_files ]
				}.transpose()

		} else reads_for_alignment = FASTP_TRIM.out.reads

	//    FASTP.out.reads.view()
	//    FASTP.out.meta.view()
	//	reads_for_alignment.view()
	    ch_multiqc_files = ch_multiqc_files.mix(FASTP_TRIM.out.json.collect{ meta, json -> json })
	    ch_multiqc_files = ch_multiqc_files.mix(FASTP_TRIM.out.html.collect{ meta, html -> html })
	    ch_versions = ch_versions.mix(FASTP_TRIM.out.versions.first())

	} else {

		reads_for_alignment = ch_samplesheet 
			
	}


	//	reads_for_alignment.view()

	
    // STEP 1: MAPPING READS TO REFERENCE GENOME
    // First, we must calculate number of lanes for each sample (meta.n_fastq)
    // This is needed to group reads from the same sample together using groupKey to avoid stalling the workflow
    // when reads from different samples are mixed together
    reads_for_alignment.map { meta, reads ->
            [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ]
        }
        .groupTuple()
        .map { meta, reads ->
            meta + [ n_fastq: reads.size() ] // We can drop the FASTQ files now that we know how many there are
        }
        .set { reads_grouping_key } 

	seq_center = "DNALINK"
	seq_platform = "ILLUMINA"
	fasta_file= params.fasta

	reads_for_alignment.map { meta, reads ->
	
		def CN         = seq_center ? "CN:${seq_center}\\t" : ''
		def flowcell   = flowcellLaneFromFastq( reads[0] )
		def read_group = "\"@RG\\tID:${flowcell}.${meta.id}\\t${CN}PU:lane1\\tSM:${meta.id}\\tLB:${meta.id}\\tDS:${fasta_file}\\tPL:${seq_platform}\""
		
		[ meta + [ read_group: read_group.toString() ], reads ]
	}.set{ reads_for_alignment }


	reads_for_alignment.view()
	// reads_grouping_key.view() //  [n_fastq:6]

    // reads will be sorted
    sort_bam = true
    FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON(reads_for_alignment, index_alignment, sort_bam, fasta, fasta_fai)

	ch_versions = ch_versions.mix(FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.versions)

	FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bam.view()


	// STEP : MARKDUP.

    // Grouping the bams from the same samples not to stall the workflow
    // Use groupKey to make sure that the correct group can advance as soon as it is complete
    // and not stall the workflow until all reads from all channels are mapped
    bam_mapped = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bam
        .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
        .filter { meta1, bam, meta2 -> meta1.sample == meta2.sample }
        // Add n_fastq and other variables to meta
        .map { meta1, bam, meta2 ->
            [ meta1 + meta2, bam ]
        }
        // Manipulate meta map to remove old fields and add new ones
        .map { meta, bam ->
            [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size') + [ data_type: 'bam', id: meta.id ], bam ]
        }
        // Create groupKey from meta map
        .map { meta, bam ->
            [ groupKey( meta, meta.n_fastq), bam ]
        }
        // Group
        .groupTuple()

    bai_mapped = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bai
        .combine(reads_grouping_key) // Creates a tuple of [ meta, bai, reads_grouping_key ]
        .filter { meta1, bai, meta2 -> meta1.sample == meta2.sample }
        // Add n_fastq and other variables to meta
        .map { meta1, bai, meta2 ->
            [ meta1 + meta2, bai ]
        }
        // Manipulate meta map to remove old fields and add new ones
        .map { meta, bai ->
            [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size') + [ data_type: 'bai', id: meta.sample ], bai ]
        }
        // Create groupKey from meta map
        .map { meta, bai ->
            [ groupKey( meta, meta.n_fastq), bai ]
        }
        // Group
        .groupTuple()

/*
	// bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
	BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)
	BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)

	// Gather used softwares versions
	ch_versions = ch_versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
	ch_versions = ch_versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
*/




	cram_for_markduplicates = bam_mapped
	intervals_for_preprocessing = Channel.value([ [ id:'null' ], [] ])
	
	
	cram_for_markduplicates.view()

    BAM_MARKDUPLICATES(
		cram_for_markduplicates,
		fasta,
		fasta_fai,
		intervals_for_preprocessing)

	cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

	
	// Gather QC reportsm_for_markduplicates
	ch_multiqc_filess = ch_multiqc_files.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> [ report ] })
	// Gather used softwares versions
	ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)


	// ch_md_cram_for_restart contains either:
	// - crams from markduplicates
	// - crams from sentieon_dedup
	// - crams from markduplicates_spark
	// - crams from input step markduplicates --> from the converted ones only?
	ch_md_cram_for_restart = Channel.empty().mix(cram_markduplicates_no_spark)
	    // Make sure correct data types are carried through
	    .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }


	ch_md_cram_for_restart.view()

	
	// If params.save_output_as_bam, then convert CRAM files to BAM
	CRAM_TO_BAM(ch_md_cram_for_restart, fasta, fasta_fai)
	ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)


	CRAM_TO_BAM.out.cram.view()




    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: false
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}



