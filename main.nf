#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020

This software is a computer program whose purpose is to
analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms
of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful,
but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards
their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge
of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                               Variant calling pipeline
========================================================================================
#### Homepage / Documentation
https://gitlab.curie.fr/data-analysis/tumospec/-/tree/tumospec
----------------------------------------------------------------------------------------
*/

// File with text to display when a developement version is used
devMessageFile = file("$baseDir/assets/devMessage.txt")

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
     log.info devMessageFile.text
  }

  log.info """
  v${workflow.manifest.version}
  ======================================================================

  Usage:
  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile conda
  nextflow run main.nf --samplePlan samplePlan --genome 'hg19' -profile conda

  Mandatory arguments:
    --reads [file]                Path to input data (must be surrounded with quotes)
    --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
    --genome [str]                Name of genome reference
    -profile [str]                Configuration profile to use. test / conda / multiconda / path / multipath / singularity / docker / cluster (see below)
    --snpeffDb [str]              Directory to the snpEff databases

  Inputs:
    --design [file]               Path to design file for extended analysis  
    --singleEnd [bool]            Specifies that the input is single-end reads

  Skip options: All are false by default
    --skipSoftVersion [bool]      Do not report software version
    --skipMultiQC [bool]          Skips MultiQC

  Other options:
    --outDir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
 
  =======================================================
  Available Profiles

    -profile test                Set up the test dataset
    -profile conda               Build a single conda for with all tools used by the different processes before running the pipeline
    -profile multiconda          Build a new conda environment for each tools used by the different processes before running the pipeline
    -profile path                Use the path defined in the configuration for all tools
    -profile multipath           Use the paths defined in the configuration for each tool
    -profile docker              Use the Docker containers for each process
    -profile singularity         Use the singularity images for each process
    -profile cluster             Run the workflow on the cluster, instead of locally

  """.stripIndent()
}


/***********************************************************************************************************
                                    * SET UP CONFIGURATION VARIABLES *
***********************************************************************************************************/

nextflow.enable.dsl=1

// Show help message
if (params.help){
  helpMessage()
  exit 0
}

// Configurable reference genomes

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}


params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwaIndex = params.genome ? params.genomes[ params.genome ].bwaIndex ?: false : false
params.geneBed = params.genome ? params.genomes[ params.genome ].geneBed ?: false : false
params.nameSnpeff = params.genome ? params.genomes[ params.genome ].nameSnpeff ?: false : false
params.gnomad = params.genome ? params.genomes[ params.genome ].gnomad ?: false : false
params.dbNSFP = params.genome ? params.genomes[ params.genome ].dbNSFP ?: false : false


GenePositionCh = params.genes_tumospec  ? Channel.fromPath(params.genes_tumospec, checkIfExists: true).collect() : Channel.empty()
VariantTypeCh = params.var_types  ? Channel.fromPath(params.var_types, checkIfExists: true).collect() : Channel.empty()
SampleNameCh = params.samples_plan  ? Channel.fromPath(params.samples_plan, checkIfExists: true).collect() : Channel.empty()
ExonCh = params.exons_tumospec  ? Channel.fromPath(params.exons_tumospec, checkIfExists: true).collect() : Channel.empty()
caddConfCh = params.cadd_conf  ? Channel.fromPath(params.cadd_conf, checkIfExists: true).collect() : Channel.empty()



if ( params.fasta ){
Channel.fromPath(params.fasta)
  .ifEmpty { exit 1, "Reference annotation not found: ${params.fasta}" }
  .set { fastaCh }
}else{
  fastaCh = Channel.empty()
}


// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

// Stage config files
MultiqcConfigCh = Channel.fromPath(params.multiqcConfig)
outputDocsCh = file("$baseDir/docs/output.md", checkIfExists: true)
outputDocsImagesCh = file("$baseDir/docs/images/", checkIfExists: true)


/***********************************************************************************************************
                                             * CHANNELS *
***********************************************************************************************************/

// Validate inputs
if ((params.reads && params.samplePlan) || (params.readPaths && params.samplePlan)){
  exit 1, "Input reads must be defined using either '--reads' or '--samplePlan' parameters. Please choose one way."
}

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { metadataCh }
}else{
  metadataCh = Channel.empty()
}

// Create a channel for input read files
if(params.samplePlan){
  if(params.singleEnd){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .set { rawReadsCh }
  }else{
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .set { rawReadsCh }
   }
  params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .set { rawReadsCh }
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .set { rawReadsCh }
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { rawReadsCh }
}

// Make sample plan if not available
if (params.samplePlan){
  Channel
    .fromPath(params.samplePlan)
    .into{ samplePlanCh; samplePlanCheckCh }
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
      .from(params.readPaths)
      .collectFile() {
        item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .into{ samplePlanCh; samplePlanCheckCh }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .into{ samplePlanCh; samplePlanCheckCh }
  }
}else{
  if (params.singleEnd){
    Channel
      .fromFilePairs( params.reads, size: 1 )
      .collectFile() {
        item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .into { samplePlanCh; samplePlanCheckCh }
  }else{
    Channel
      .fromFilePairs( params.reads, size: 2 )
      .collectFile() {
        item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
      }
      .into { samplePlanCh; samplePlanCheckCh }
   }
}


/***********************************************************************************************************
                                          * Design file *
***********************************************************************************************************/

// TODO - UPDATE BASED ON THE HEADER OF YOUR DESIGN

if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .into { designCheckCh ; designCh }

  designCh
    .splitCsv(header:true)
    .map { row ->
      return [ row.SAMPLEID, row.CONTROLID, row.SAMPLENAME, row.GROUP, row.PEAKTYPE ]
     }
    .set { designCh }
}else{
  designCheckCh = Channel.empty()
  designCh = Channel.empty()
}


/***********************************************************************************************************
                                         * Header log info *
***********************************************************************************************************/

if ("${workflow.manifest.version}" =~ /dev/ ){
   log.info devMessageFile.text
}

log.info """=======================================================

workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]

summary['Pipeline Name']  = 'Variant calling pipeline'
summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Container Engine'] = workflow.containerEngine
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/***********************************************************************************************************
                                           * Processes *
***********************************************************************************************************/


/***********************************************************************************************************
                                     * Control quality with FastQC *
***********************************************************************************************************/

// To split a channel in two channels
rawReadsCh.into { rawReadsFastqcCh; rawReadsTrimGaloreCh }


process fastqc {
  label 'fastqc'
  label 'lowMem'
  label 'lowCpu'

  tag "${sample_id}"
  publishDir "${params.outDir}/01_Fastqc", mode: 'copy'

  input:
  set val(sample_id), file(reads) from rawReadsFastqcCh

  output:
  file("*_fastqc.{zip,html}") into FastqcReportCh
  file "v_fastqc.txt" into FastqcVersionCh
  
  script:
  """
  fastqc -q --threads ${task.cpus} ${reads}
  fastqc --version &> v_fastqc.txt 2>&1 || true
  """
}


/**********************************************************************************************************
                                      * Trimming with TrimGalore *
***********************************************************************************************************/
// Trimming of the adapters sequences and the low quality ends                                                                                           
// -q Trim low-quality ends from reads in addition to adapter removal (by default phred score at 20)                                                     
// --length Discard reads that became shorter than length INT because of either quality or adapter trimming (default 20 bp) 

process trimgalore {
  label 'trimgalore'
  label 'lowMem'
  label 'highCpu'

  publishDir "${params.outDir}/02_TrimGalore", mode: 'copy'

  input:
  set val(sample_id), file(reads) from rawReadsTrimGaloreCh

  output:
  set val(sample_id), file("*fastq.gz") into trimReadsTrimGaloreCh
  file("*_trimming_report.txt") into TrimGaloreReportCh
  file "v_trimgalore.txt" into TrimGaloreVersionCh
  
  script:
  """
  trim_galore -q ${params.baseQualityScoreEnd} --length ${params.fragmentSize} --paired --gzip ${reads} --basename ${sample_id} --cores ${task.cpus} 
  mv ${sample_id}_val_1.fq.gz ${sample_id}_R1_trimmed.fastq.gz
  mv ${sample_id}_val_2.fq.gz ${sample_id}_R2_trimmed.fastq.gz
  trim_galore --version &> v_trimgalore.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                   * Control quality with FastQC after trimming *
***********************************************************************************************************/

// To split a channel in two channels
trimReadsTrimGaloreCh.into { trimReadsFastqcCh; trimReadsBwaCh }


process fastqc_trimmed {
  label 'fastqc'
  label 'lowMem'
  label 'lowCpu'

  publishDir "${params.outDir}/03_FastqcTrimmed", mode: 'copy'

  input:
  set val(sample_id), file(reads) from trimReadsFastqcCh

  output:
  file("*_fastqc.{zip,html}") into trimFastqcReportCh
  
  script:
  """
  fastqc -q --threads ${task.cpus} $reads
  """
}


/***********************************************************************************************************
                         * Alignment with BWA MEM and produce sort BAM with Samtools *
***********************************************************************************************************/

process bwa_mem {
  label 'bwa'
  label 'extraMem'
  label 'extraCpu'

  publishDir "${params.outDir}/04_Bam", mode: 'copy'

  input:
  set val(sample_id), file(reads) from trimReadsBwaCh

  output:
  set val(sample_id), file("*.sort.bam") into SortBamCh
  set val(sample_id), file("*.sort.bam.bai") into IndexSortBamCh
  file "v_bwa.txt" into BwaVersionCh


  script:
  readGroup = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:LB_${sample_id}\\tPL:Illumina" 
  """
  bwa mem -t ${task.cpus} \
  -R \"${readGroup}\" \
  ${params.bwaIndex} ${reads} \
  | samtools view -@ ${task.cpus} -bh \
  | samtools sort -@ ${task.cpus} > ${sample_id}.sort.bam
  samtools index ${sample_id}.sort.bam
  bwa &> v_bwa.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                      * Filtering alignment with Samtools *
***********************************************************************************************************/
// Filtering alignment file on mapping quality and delete non aligned reads                                                                           
// -b output in bam format                                                                                                                            
// -h include header in output                                                                                                                        
// -q skip alignment with MAPQ smaller than a value                                                                                                   
// -F filter by flag. flag=4 get only the mapped data 


process samtools_filtering {
  label 'samtools'
  label 'lowMem'
  label 'lowCpu'

  publishDir "${params.outDir}/04_Bam", mode: 'copy'

  input:
  set val(sample_id), file(bam) from SortBamCh

  output:
  set val(sample_id), file("*.sort.filt.bam") into FiltBamCh
  set val(sample_id), file("*.sort.filt.bam.bai") into IndexFiltBamCh
  file('v_samtools.txt') into SamtoolsVersionCh

  script:
  """
  samtools view -bh -F ${params.flag} -q ${params.mappingQuality} ${bam} > ${sample_id}.sort.filt.bam
  samtools index ${sample_id}.sort.filt.bam  
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                   * Intersection with the BED file with Bedtools *
***********************************************************************************************************/

process bedtools_intersect {
  label 'bedtools'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/04_Bam", mode: 'copy'

  input:
  set val(sample_id), file(bam) from FiltBamCh

  output:
  set val(sample_id), file("*.sort.filt.ontarget.bam") into TargetBamCh
  set val(sample_id), file("*.sort.filt.ontarget.bam.bai") into IndexTargetBamCh
  file('v_bedtools.txt') into BedtoolsVersionCh

  script:
  """
  bedtools intersect -a ${bam} -b ${params.geneBed} > ${sample_id}.sort.filt.ontarget.bam 
  samtools index ${sample_id}.sort.filt.ontarget.bam
  bedtools --version &> v_bedtools.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                   * Coverage analysis with Samtools mpileup *
***********************************************************************************************************/
// -A Do not skip anomalous read pairs in variant calling                                                                  
// -B Disable probabilistic realignment for the computation of base alignment quality (BAQ)                                
// -d At a position, read maximally INT reads per input file (min value : 8000)                                            
// -q Minimum mapping quality for an alignment to be used [0]                                                              
// -Q Minimum base quality for a base to be considered [13]                                                                
// -x Disable read-pair overlap detection            

// To split a channel in two channels
TargetBamCh.into { TargetBamMpileupCh; TargetBamGatkCh }

process samtools_mpileup {
  label 'samtools'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/05_Mpileup", mode: 'copy'

  input:
  set val(sample_id), file(bam) from TargetBamMpileupCh

  output:
  set val (sample_id), file("*.mpileup") into MpileupReportCh
  set val (sample_id), file("*.mpileup.csv") into MpileupCsvCh

  script:
  """
  samtools mpileup -A -B -x \
  -d ${params.maxReads} \
  -q ${params.mappingQuality} \
  -l ${params.geneBed} \
  -Q ${params.baseQualityScore} \
  -f ${params.fasta} ${bam} > ${sample_id}.mpileup

  awk -F " " '{print \$1 "," \$2 "," \$3 "," \$4}' ${sample_id}.mpileup > ${sample_id}.mpileup.csv
  
  """
}

 
/***********************************************************************************************************
                                         * Coverage analysis with Python *
***********************************************************************************************************/


process python_coverage {
  label 'python'
  label 'lowMem'
  label 'medCpu'

publishDir "${params.outDir}/05_Mpileup", mode: 'copy'

  input:
  set val(sample_id), file(mpileup) from MpileupCsvCh
  file(gene_position) from GenePositionCh.collect()
  
  output:
  set val(sample_id), file("*.coverage.mpileup.python.csv") into PythonMpileupCsvCh

  script:
  """
  coverage_mpileup.py ${mpileup} ${gene_position}
  """
}



/***********************************************************************************************************
                                          * Variant calling with GATK *
***********************************************************************************************************/
// -I : input BAM/SAM/CRAM file containing reads                                                                           
// -L : One or more genomic intervals over which to operate                                                                
// -R : Reference sequence file                                                                                  
// --min-base-quality-score or -mbq : Minimum base quality required to consider a base for calling. Default value : 10     
// --minimum-mapping-quality : Minimum mapping quality required to consider a base for calling. Default value : 10         
//  -stand-call-conf : The minimum phred-scaled confidence threshold at which variants should be called. Default value : 30
// --emit-ref-confidence or -ERC : Mode for emitting reference confidence scores                                           
// -O : output File to which variants should be written 


process gatk_variant_calling {
  label 'gatk'
  label 'medMem'
  label 'highCpu'

publishDir "${params.outDir}/06_GatkVcf", mode: 'copy'

  input:
  set val(sample_id), file(bam) from TargetBamGatkCh

  output:
  set val(sample_id), file("*.GATK.vcf") into VcfCh
  set val(sample_id), file("*.GATK.vcf.idx") into IndexVcfCh
  file('v_gatk.txt') into GatkVersionCh

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}G" \
  HaplotypeCaller \
  -I ${bam} \
  -L ${params.geneBed} \
  -R ${params.fasta} \
  --min-base-quality-score ${params.baseQualityScore} \
  --minimum-mapping-quality ${params.mappingQuality} \
  -stand-call-conf ${params.variantQuality} \
  -ERC NONE \
  -O ${sample_id}.GATK.vcf
  gatk HaplotypeCaller --help &> v_gatk.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                        * Variant filtering with GATK *
***********************************************************************************************************/
// -R : reference Reference sequence                                                                                                                  
// -V : variant a VCF file containing variants                                                                                                        
// --filter-name : Names to use for the list of filters                                                                                               
// --filter-expression : One or more expressions used with INFO fields to filter

process gatk_variant_filtering {
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/06_GatkVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from VcfCh

  output:
  set val(sample_id), file("*.GATK.prefilt.vcf") into PrefiltVcfCh
  set val(sample_id), file("*.GATK.prefilt.vcf.idx") into IndexPrefiltVcfCh

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}G" \
  VariantFiltration \
  -R ${params.fasta} \
  -V ${vcf} \
  -O ${sample_id}.GATK.prefilt.vcf \
  --filter-name "${params.qd_filter}" \
  -filter "${params.qd}" \
  --filter-name "${params.fs_filter}" \
  -filter "${params.fs}" \
  --filter-name "${params.sor_filter}" \
  -filter "${params.sor}" \
  --filter-name "${params.mq_filter}" \
  -filter "${params.mq}" \
  --filter-name "${params.mqrs_filter}" \
  -filter "${params.mqrs}" \
  --filter-name "${params.rprs_filter}" \
  -filter "${params.rprs}"
  """
}


/***********************************************************************************************************
                                      * Select variants after filtering with GATK *
***********************************************************************************************************/

process gatk_variant_selection {
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/06_GatkVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from PrefiltVcfCh

  output:
  set val(sample_id), file("*.GATK.filt.vcf") into FiltVcfCh
  set val(sample_id), file("*.GATK.filt.vcf.idx") into IndexFiltVcfCh

script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}G" \
  SelectVariants \
  -R ${params.fasta} \
  -V ${vcf} \
  --exclude-filtered \
  -O ${sample_id}.GATK.filt.vcf
  """
}


/***********************************************************************************************************
                                     * Split multi allelic sites with Bcftools *
***********************************************************************************************************/

process bcftools_split {
  label 'bcftools'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/06_GatkVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from FiltVcfCh

  output:
  set val(sample_id), file("*.GATK.filt.multi.vcf") into MultiFiltVcfCh
  file('v_bcftools.txt') into BcftoolsVersionCh

script:
  """
  bcftools norm -m - \
  -f ${params.fasta} \
  -O v -o ${sample_id}.GATK.filt.multi.vcf \
  ${vcf}
  bcftools --version &> v_bcftools.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                       * Annotate with ref genome with SnpEff *
***********************************************************************************************************/

process snpeff_annotate {
  label 'snpeff'
  label 'medMem'
  label 'lowCpu'

publishDir "${params.outDir}/07_SnpEffVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from MultiFiltVcfCh

  output:
  set val(sample_id), file("*.SNPEFF.ann.vcf") into AnnVcfCh
  file("*.snpeff_report.html") into SnpeffReportCh
  file("*.snpEff.summary.csv") into SnpeffReportCsvCh
  file('v_snpeff.txt') into SnpeffVersionCh

script:
  """
  snpEff -Xmx${task.memory.toGiga()}g \
  -v ${params.nameSnpeff} \
  -s ${sample_id}.snpeff_report.html \
  -csvStats ${sample_id}.snpEff.summary.csv \
  -dataDir ${params.snpeffDb} \
  ${vcf} \
  > ${sample_id}.SNPEFF.ann.vcf
  snpEff -help &> v_snpeff.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                    * Annotate with Gnomad database with SnpSift *
***********************************************************************************************************/

process snpsift_annotate {
  label 'snpsift'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/07_SnpEffVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from AnnVcfCh

  output:
  set val(sample_id), file("*.SNPEFF.ann.gnomad.vcf") into GnomadAnnVcfCh
  file('v_snpsift.txt') into SnpSiftVersionCh

script:
  """
  SnpSift -Xmx${task.memory.toGiga()}g annotate \
  -info ${params.ethnies} \
  ${params.gnomad} \
  ${vcf} \
  > ${sample_id}.SNPEFF.ann.gnomad.vcf
  SnpSift annotate -h &> v_snpsift.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                    * Annotate with CADD score with VcfAnno *
***********************************************************************************************************/


process cadd_annotate {
  label 'vcfanno'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/07_SnpEffVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from GnomadAnnVcfCh
  file(caddConf) from caddConfCh

  output:
  set val(sample_id), file("*.SNPEFF.ann.gnomad.cadd.vcf") into CaddGnomadAnnVcfCh
  file('v_vcfanno.txt') into VcfAnnoVersionCh

script:
  """
  vcfanno -p ${task.cpus} \
  ${caddConf} \
  ${vcf} \
  > ${sample_id}.SNPEFF.ann.gnomad.cadd.vcf
  vcfanno -help &> v_vcfanno.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                    * Annotate with dbNSFP with SnpSift *
***********************************************************************************************************/

process snpsift_dbNSFP {
  label 'snpsift'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/07_SnpEffVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from CaddGnomadAnnVcfCh

  output:
  set val(sample_id), file("*.SNPEFF.ann.gnomad.cadd.dbNSFP.vcf") into dbNSFPCaddGnomadAnnVcfCh

script:
  """
  SnpSift -Xmx${task.memory.toGiga()}g dbnsfp \
  -db ${params.dbNSFP} \
  -f ${params.scores} \
  ${vcf} \
  > ${sample_id}.SNPEFF.ann.gnomad.cadd.dbNSFP.vcf
  """
}


/***********************************************************************************************************
                                            * Compress VCF file *
***********************************************************************************************************/

// To split a channel in two channels
dbNSFPCaddGnomadAnnVcfCh.into { TobgzipdbNSFPCaddGnomadAnnVcfCh; TosnpsiftdbNSFPCaddGnomadAnnVcfCh }

process compress_vcf {
  label 'bgzip'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/07_SnpEffVcf", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from TobgzipdbNSFPCaddGnomadAnnVcfCh

  output:
  set val(sample_id), file("*.SNPEFF.ann.gnomad.cadd.dbNSFP.vcf.gz") into CompressVcfCh
  set val(sample_id), file("*.SNPEFF.ann.gnomad.cadd.dbNSFP.vcf.gz.tbi") into IndexCompressVcfCh

script:
  """
  bgzip -c ${vcf} \
  > ${sample_id}.SNPEFF.ann.gnomad.cadd.dbNSFP.vcf.gz

  tabix -p vcf ${sample_id}.SNPEFF.ann.gnomad.cadd.dbNSFP.vcf.gz
  """
}


/***********************************************************************************************************
                                        * Convert VCF into TXT file with SnpSift *
***********************************************************************************************************/

process snpsift_convert {
  label 'snpsift'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/08_Txt", mode: 'copy'

  input:
  set val(sample_id), file(vcf) from TosnpsiftdbNSFPCaddGnomadAnnVcfCh

  output:
  file("*.SNPEFF.ann.gnomad.cadd.dbNSFP.recode.txt") into TxtVcfCh

script:
  """
   cat ${vcf} \
   | vcfEffOnePerLine.pl \
   | SnpSift -Xmx${task.memory.toGiga()}g \
   extractFields -s "," -e "." - \
   ${params.col} \
   > ${sample_id}.SNPEFF.ann.gnomad.cadd.dbNSFP.recode.txt
  """
}

// To select NM add this line before SnpSift command and after vcfEffonePerLine
//  | grep -E ${params.nm} \


/***********************************************************************************************************
                                          * Filtering variants with Python *
***********************************************************************************************************/

process python_filtering_variants {
  label 'python'
  label 'lowMem'
  label 'lowCpu'

publishDir "${params.outDir}/09_FilteringVariantsPyhton", mode: 'copy'

  input:
  file(txt) from TxtVcfCh.collect()
  file(variants_type) from VariantTypeCh
  file(sample_name) from SampleNameCh
  file(exons) from ExonCh

  output:
  file("df_variants.csv") into VariantsCh
  file("df_variants_formatted.csv") into FormattedVariantsCh
  file("number_of_variants.txt") into NumberVariantsCh

script:
  """
  filtering_variants.py ${txt} ${variants_type} ${sample_name} ${exons}
  """
}


/***********************************************************************************************************
                                                 * MultiQC *
***********************************************************************************************************/

process getSoftwareVersions{
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/10_SoftwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file 'v_fastqc.txt' from FastqcVersionCh.first().ifEmpty([])
  file 'v_trimgalore.txt' from TrimGaloreVersionCh.first().ifEmpty([])
  file 'v_bwa.txt' from BwaVersionCh.first().ifEmpty([])
  file 'v_samtools.txt' from SamtoolsVersionCh.first().ifEmpty([])
  file 'v_bedtools.txt' from BedtoolsVersionCh.first().ifEmpty([])
  file 'v_gatk.txt' from GatkVersionCh.first().ifEmpty([])
  file 'v_bcftools.txt' from BcftoolsVersionCh.first().ifEmpty([])
  file 'v_snpeff.txt' from SnpeffVersionCh.first().ifEmpty([])
  file 'v_snpsift.txt' from SnpSiftVersionCh.first().ifEmpty([])
  file 'v_vcfanno.txt' from VcfAnnoVersionCh.first().ifEmpty([])

  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYamlCh

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}


process workflowSummaryMqc {
  when:
  !params.skipMultiQC

  output:
  file 'workflow_summary_mqc.yaml' into workflowSummaryYamlCh

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: "${workflow.manifest.homePage}"
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}


process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'medMem'

  publishDir "${params.outDir}/11_MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from samplePlanCh.collect()
  file MultiqcConfig from MultiqcConfigCh
  file ('fastqc/*') from FastqcReportCh.collect().ifEmpty([])
  file ('trimming/*') from TrimGaloreReportCh.collect().ifEmpty([])
  file ('fastqc/*') from trimFastqcReportCh.collect().ifEmpty([])
  file ('snpEff/*') from SnpeffReportCsvCh.collect().ifEmpty([])
  file metadata from metadataCh.ifEmpty([])
  file ('softwareVersions/*') from softwareVersionsYamlCh.collect().ifEmpty([])
  file ('workflowSummary/*') from workflowSummaryYamlCh.collect()

  output: 
  file splan
  file "*report.html" into MultiqcReportCh
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  modulesList = "-m custom_content -m fastqc -m cutadapt -m snpeff"
  """
  mqc_header.py --splan ${splan} --name "${workflow.manifest.name}" --version "${workflow.manifest.version}" ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $MultiqcConfig $modulesList
  multiqc --help &> v_multiqc.txt 2>&1 || true
  """
}


/***********************************************************************************************************
                                             * Sub-routines *
***********************************************************************************************************/

process checkDesign{
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.SummaryDir}/", mode: 'copy'

  when:
  params.design

  input:
  file design from designCheckCh
  file samplePlan from samplePlanCheckCh

  script:
  optSE = params.singleEnd ? "--singleEnd" : ""
  """
  checkDesign.py -d $design -s $samplePlan ${optSE}
  """
}


process outputDocumentation {
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.summaryDir}/", mode: 'copy'

  input:
  file outputDocs from outputDocsCh
  file images from outputDocsImagesCh

  output:
  file "resultsDescription.html"

  script:
  """
  markdown_to_html.py $outputDocs -o resultsDescription.html
  """
}

workflow.onComplete {

  // pipelineReport.html
  def reportFields = [:]
  reportFields['pipeline'] = workflow.manifest.name
  reportFields['version'] = workflow.manifest.version
  reportFields['runName'] = customRunName ?: workflow.runName
  reportFields['success'] = workflow.success
  reportFields['dateComplete'] = workflow.complete
  reportFields['duration'] = workflow.duration
  reportFields['exitStatus'] = workflow.exitStatus
  reportFields['errorMessage'] = (workflow.errorMessage ?: 'None')
  reportFields['errorReport'] = (workflow.errorReport ?: 'None')
  reportFields['commandLine'] = workflow.commandLine
  reportFields['projectDir'] = workflow.projectDir
  reportFields['summary'] = summary
  reportFields['summary']['Date Started'] = workflow.start
  reportFields['summary']['Date Completed'] = workflow.complete
  reportFields['summary']['Pipeline script file path'] = workflow.scriptFile
  reportFields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) reportFields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) reportFields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) reportFields['summary']['Pipeline Git branch/tag'] = workflow.revision

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/workflowOnCompleteTemplate.txt")
  def txtTemplate = engine.createTemplate(tf).make(reportFields)
  def reportTxt = txtTemplate.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/workflowOnCompleteTemplate.html")
  def htmlTemplate = engine.createTemplate(hf).make(reportFields)
  def reportHtml = htmlTemplate.toString()

  // Write summary HTML to a file
  def outputSummaryDir = new File( "${params.summaryDir}/" )
  if( !outputSummaryDir.exists() ) {
    outputSummaryDir.mkdirs()
  }
  def outputHtmlFile = new File( outputSummaryDir, "pipelineReport.html" )
  outputHtmlFile.withWriter { w -> w << reportHtml }
  def outputTxtFile = new File( outputSummaryDir, "pipelineReport.txt" )
  outputTxtFile.withWriter { w -> w << reportTxt }

  // workflowOnComplete file
  File woc = new File("${params.outDir}/workflowOnComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'
  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  // final logs
  if(workflow.success){
    log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}
