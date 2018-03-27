#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         drax
========================================================================================
 drax Analysis Pipeline. Started 2018-03-05.
 #### Homepage / Documentation
 https://github.com/will-rowe/drax
 #### Authors
 Will Rowe will-rowe <will.rowe@stfc.ac.uk> - https://will-rowe.github.io>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     drax v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    drax --reads '*_R{1,2}.fastq.gz' --refData ./DRAX-files

    Mandatory arguments:
      --reads           Path to input data (must be surrounded with quotes)
      --refData         Path to reference data (run `drax get` to download this)

    QC options:
      --singleEnd       Specifies that the input is single end reads
      --qual            Quality cut-off to use in QC workflow (deduplication, trimming etc.)
      --decontaminate   Peform read subtraction using masked HG19 reference (default = true)

    Other options:
      --outdir          The output directory where the results will be saved
      --email           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      -profile          Specify docker or standard (default)
      --max_cpus        Number of CPUs to use (default = max available)
      --max_memory      Amount of RAM to use (default = 32 GB)
    """.stripIndent()
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    CONFIG
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1'

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.singleEnd = false
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.reads = ""
params.refData  = ""
params.outdir = './drax-results'
params.email = false
params.plaintext_email = false
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
params.qual = 20
params.decontaminate = true

// Validate inputs
if ( params.reads  == "" ) exit 1, "Must provide reads ('some/path/*_R{1,2}.fastq.gz')"
if ( params.refData  == "" ) exit 1, "Must provide reference data (--refData), use `drax get` to download it"
refData = file(params.refData)
if ( !refData.exists() ) exit 1, "Supplied refData  does not exist!"
if ( !refData.isDirectory() ) exit 1, "Supplied refData is not a directory!"
// TODO: check here that the refData dir has the correct contents
refDataLink = file("${workflow.workDir}/ref-data")
if( !refDataLink.exists() ) {
    refData.mklink(refDataLink)
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info "========================================="
log.info " drax v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Output Dir']   = params.outdir
summary['Working Dir']  = workflow.workDir
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
summary['Reference Data'] =   params.refData
summary['Quality Cut-off']  =   params.qual
summary['Decontaminate'] = params.decontaminate
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// create some additional log files
logDir = file("${params.outdir}/logs")
if( !logDir.exists() ) {
    logCheck = logDir.mkdir()
    if ( !logCheck ) exit 1, "Cannot create log directory: $logDir"
}
logFileForDeduplicate = file(logDir + "/deduplicate.log")
logFileForTrimming = file(logDir + "/trimming.log")
logFileForReadSubtraction = file(logDir + "/readSubtraction.log")


/*
 * Parse software version numbers
 */
process get_software_versions {
    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * Create channels
 */
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { input_data }

input_data.into { read_files_fastqc; read_files_to_deduplicate }

cpus = params.max_cpus

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    QUALITY  CONTROL
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc-initial-check", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads -t $cpus
    """
}


/*
 * Deduplicate
 */
process deduplicate {
    tag "$sampleID"

    input:
    set sampleID, file(reads) from read_files_to_deduplicate

	output:
    set sampleID, file("${sampleID}*.deduplicated.fq.gz") into read_files_to_trim
    file("${sampleID}.deduplicate.log") into logDeduplicate

    script:
	"""
    # log some stuff
    echo "------------------------------------------------------" >> ${sampleID}.deduplicate.log
    echo "SAMPLE: ${sampleID}" >> ${sampleID}.deduplicate.log
    echo "------------------------------------------------------" >> ${sampleID}.deduplicate.log

    # set up the command
    if [ \"$params.singleEnd\" = \"false\" ]; then
		dedupeCMD=\"clumpify.sh dedupe in1=${reads[0]} in2=${reads[1]} out1=${sampleID}_R1.deduplicated.fq.gz out2=${sampleID}_R2.deduplicated.fq.gz subs=0 threads=${cpus}\"
    else
        dedupeCMD=\"clumpify.sh dedupe in1=${reads[0]} out1=${sampleID}_R1.deduplicated.fq.gz subs=0 qin=${params.qual} threads=${cpus}\"
	fi

    # run the command
    \$dedupeCMD 2>&1 | tee .tmp

    # parse the command output and log more stuff
    duplicatesFound=\$(grep \"Duplicates Found:\" .tmp | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
    readsIn=\$(grep \"Reads In:\" .tmp | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
    remainingReads=\$((\$readsIn-\$duplicatesFound))
    percentage=\$(echo \$remainingReads \$readsIn | awk '{print \$1/\$2*100}' )
    sed -n '/Reads In:/,/Total time:/p' .tmp >> ${sampleID}.deduplicate.log
    printf "\n\$percentage%% of reads retained.\n\n" >> ${sampleID}.deduplicate.log
	"""
}


/*
 * Trimming
 */
process trimming {
    tag "$sampleID"

    input:
    set sampleID, file(reads) from read_files_to_trim

	output:
    file("${sampleID}.trimming.log") into logTrimming
    set sampleID, file("${sampleID}*.deduplicated.trimmed.fq.gz") into read_files_to_readSubtraction

    script:
	"""
    # log some stuff
    echo "------------------------------------------------------" >> ${sampleID}.trimming.log
    echo "SAMPLE: ${sampleID}" >> ${sampleID}.trimming.log
    echo "------------------------------------------------------" >> ${sampleID}.trimming.log

    # set up the command
    if [ \"$params.singleEnd\" = \"false\" ]; then
		fastpCMD=\"fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_R1.deduplicated.trimmed.fq.gz -O ${sampleID}_R2.deduplicated.trimmed.fq.gz -M ${params.qual} --cut_by_quality5 -w ${cpus}\"
    else
        fastpCMD=\"fastp -i ${reads[0]} -o ${sampleID}_R1.deduplicated.trimmed.fq.gz -M ${params.qual} --cut_by_quality5 -w ${cpus}\"
	fi

    # run the command
    \$fastpCMD 2>&1 | tee .tmp

    # parse the command output and log more stuff
    sed -n '/Read1 before filtering/,/bases trimmed due to adapters/p' .tmp >> ${sampleID}.trimming.log
	"""
}


/*
 * ReadSubtraction
 */
process readSubtraction {
    tag "$sampleID"
    publishDir "${params.outdir}/clean_data", mode: 'copy', pattern: "*clean*"

    input:
    set sampleID, file(reads) from read_files_to_readSubtraction

	output:
    file("${sampleID}.readSubtraction.log") into logReadSubtraction
    file("*_clean.fq.gz")
    file("${sampleID}_clean.stats") into quality_filtered_stats
    set sampleID, file("${sampleID}*_clean.fq.gz") into quality_filtered_reads
    set sampleID, file("${sampleID}*_clean.fq.gz") into quality_filtered_reads_copy
    set sampleID, file("${sampleID}*_clean.fq.gz") into quality_filtered_reads_for_kaiju

    script:
	"""
    # log some stuff
    echo "------------------------------------------------------" >> ${sampleID}.readSubtraction.log
    echo "SAMPLE: ${sampleID}" >> ${sampleID}.readSubtraction.log
    echo "------------------------------------------------------" >> ${sampleID}.readSubtraction.log

    if [ \"${params.decontaminate}\" = \"false\" ]; then
        interleaveCMD=\"interleave-fastq.py ${reads[0]} ${reads[1]}\"
        \$interleaveCMD > ${sampleID}_clean.fq
        gzip ${sampleID}_clean.fq
    else
    	if [ \"$params.singleEnd\" = \"false\" ]; then
    		readSubtractionCMD=\"bbwrap.sh -Xmx32g mapper=bbmap quickmatch fast ow=true append=t in1=${reads[0]} in2=${reads[1]} outu=${sampleID}_clean.fq.gz outm=${sampleID}_contamination.fq minid=0.97 maxindel=5 minhits=2 threads=${cpus} path=${workflow.workDir}/ref-data\"
        else
            readSubtractionCMD=\"bbwrap.sh -Xmx32g mapper=bbmap quickmatch fast ow=true append=t in1=${reads[0]} outu=${sampleID}_clean.fq.gz outm=${sampleID}_contamination.fq minid=0.97 maxindel=5 minhits=2 threads=${cpus} path=${workflow.workDir}/ref-data\"
    	fi
        # run the command
        \$readSubtractionCMD 2>&1 | tee .tmp

        # parse the command output and log more stuff
        sed -n '/Read 1 data:/,/Total time/p' .tmp >> ${sampleID}.readSubtraction.log
        cat .tmp >> ${sampleID}.readSubtraction.log
    fi

    # get some stats on the QC'd reads
    seqkitCMD=\"seqkit stats --quiet -T --threads ${cpus} ${sampleID}_clean.fq.gz\"
    \$seqkitCMD 2>&1 | tee .tmp

    # delete the header line from seqkit and send the file
    tail -n +2 .tmp > ${sampleID}_clean.stats
	"""
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    POST FILTERING QUALITY CHECK
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Post QC FastQC check
 */
process postQCcheck {
    tag "$sampleID"
    publishDir "${params.outdir}/fastqc-postqc-check", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set sampleID, file(reads) from quality_filtered_reads

    output:
    file "*_fastqc.{zip,html}" into fastqc_results2

    script:
    """
    fastqc -q $reads -t $cpus
    """
}


/*
 * MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/multiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc-initial-check/*') from fastqc_results.collect()
    file ('fastqc-postqc-check/*') from fastqc_results2.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    RESISTOME PROFILING
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
* Collect the stats from the QC'd reads and store in a single file - then create a GROOT index
*/
Channel
    .from('clean_data.stats').combine(quality_filtered_stats.flatMap())
    .collectFile(newLine: false, storeDir: params.outdir)
    .set { combined_stats_file }

process generate_groot_index {
    input:
    file(combined_stats) from combined_stats_file

    output:
    file "grootIndex" into groot_index

    script:
    """
    # get the average length column from the seqkit output
    cut -f 7 ${combined_stats} >> averageReadLength.txt
    # get the mean and stdev
    meanRL=\$(awk \'{ x+=\$1; next } END { if (x > 0) printf \"%3.0f\", x/NR}\' averageReadLength.txt)
    stdev=\$(awk \'{sum+=\$1; sumsq+=\$1*\$1} END {print sqrt(sumsq/NR - (sum/NR)^2)}\' averageReadLength.txt)

    # check that we can generate a GROOT index suitable for all samples
    cutoff=10
    if [ \$stdev -gt \$cutoff ]; then
        echo "The mean read length of the cleaned data is too variable to run GROOT on all samples!"; exit;
    fi

    # set up the command
    grootIndexCMD="groot index -i ${workflow.workDir}/ref-data/grootDB -o grootIndex -l \$meanRL -p ${cpus} -y groot-index.log"

    # run the commands
    \$grootIndexCMD 2>&1 | tee .tmp
    """
}


/*
* Combine index and clean reads channel, then run groot align
*/
toGROOT = quality_filtered_reads_copy.combine(groot_index)

process groot {
    publishDir "${params.outdir}/groot", mode: 'copy'

    input:
    set sampleID, file(reads), file(grootIndex) from toGROOT

    output:
    file "grootIndex"
    file "*.log"
    file "*.bam"
    file "*.groot-graphs"
    file "*.report"
    file "*.report" into groot_reports
    set sampleID, file(reads) into reads_for_metacherchant

    script:
    """
    # run the commands
    groot align -i ${grootIndex} -f ${reads} -p ${cpus} -y ${sampleID}.groot-align.log -o ${sampleID}.groot-graphs > out.bam
    groot report -i out.bam --lowCov -y ${sampleID}.groot-report.log > ${sampleID}-groot.report

    # notify if groot report is empty
    if [ ! -s ${sampleID}-groot.report ]; then
        echo "no ARGs were found by GROOT for ${sampleID}..."
    fi
    """
 }


 /*
 * Collect the groot reports and process them to get a list of ARGs found accross the samples
 */
 Channel
     .from('grootreports').combine(groot_reports.flatMap())
     .collectFile(newLine: false, storeDir: "${params.outdir}/groot")
     .set { combined_groot_reports }

process get_ARGs {
     input:
     file(groot_report) from combined_groot_reports

     output:
     file "groot-detected-args.fna" into detected_args

     script:
     """
     # terminate if no ARGs were found in any sample
     if [ ! -s ${groot_report} ]; then
        echo "No ARGs were found in any samples!"; exit;
     fi

     # download ARG-annot and index it
     wget https://github.com/will-rowe/groot/raw/master/db/full-ARG-databases/resfinder/resfinder.fna
     samtools faidx resfinder.fna

     # extract all the ARGs found by groot
     samtools faidx resfinder.fna `cut -f1 ${groot_report}` > groot-detected-args.fna

     # remove duplicates
     cat groot-detected-args.fna | seqkit rmdup -s -o groot-detected-args.fna
     """
}


/*
* Run metacherchant
*/
toMETACHERCHANT = reads_for_metacherchant.combine(detected_args)

process metacherchant {
   publishDir "${params.outdir}/metacherchant", mode: 'copy'

    input:
    set sampleID, file(reads), file(detectedARGs) from toMETACHERCHANT

    output:
    file "groot-detected-args.fna"
    file  "./${sampleID}"

    script:
    """
        metacherchant.sh --tool environment-finder \
            -k 31 \
            --coverage=2 \
            --maxradius 1000 \
            --reads ${reads} \
            --seq ${detectedARGs} \
            --output "./${sampleID}" \
            --work-dir "./${sampleID}/workDir" \
            -p ${cpus} \
            --trim
    """
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    TAXANOMIC PROFILING
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
* Collect the QC'd reads and run Kaiju
*/
process kaiju {
    publishDir "${params.outdir}/kaiju", mode: 'copy'

    input:
    set sampleID, file(reads) from quality_filtered_reads_for_kaiju

    output:
    file "*kaiju.out"
    file "*.html"

    script:
    """
    # run the commands
    kaiju -z ${cpus} -t ${workflow.workDir}/ref-data/nodes.dmp -f ${workflow.workDir}/ref-data/kaiju_db.fmi -i ${reads} -o ${sampleID}.kaiju.out
    kaiju2krona -t ${workflow.workDir}/ref-data/nodes.dmp -n ${workflow.workDir}/ref-data/names.dmp -i ${sampleID}.kaiju.out -o ${sampleID}.kaiju.out.krona
    ktImportText -o ${sampleID}.kaiju.out.html ${sampleID}.kaiju.out.krona
    """
 }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    CLEAN UP
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Output logs
 */
process output_logs {
    input:
    file(tolog1) from logDeduplicate.flatMap()
    file(tolog2) from logTrimming.flatMap()
    file(tolog3) from logReadSubtraction.flatMap()

    script:
    """
    cat $tolog1 >> $logFileForDeduplicate
    cat $tolog2 >> $logFileForTrimming
    cat $tolog3 >> $logFileForReadSubtraction
    """
}

/*
 * Output description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[drax] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[drax] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['software_versions'] = software_versions

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[drax] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[drax] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "========================================="
    log.info "[drax] Pipeline complete"
    log.info "[drax] Results are in: ${params.outdir}"
    log.info "========================================="
}
