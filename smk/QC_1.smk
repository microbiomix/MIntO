#!/usr/bin/env python

'''
Pre-processing of metaG and metaT data step - quality reads filtering

Authors: Carmen Saenz, Mani Arumugam
'''

import re
import os
import pandas
import glob

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule QC_1_base, QC_1_rpkg from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

######################################
# Get file suffixes
######################################

ilmn_suffix = ["1.fq.gz", "2.fq.gz"]
if (x := validate_optional_key(config, 'ILLUMINA_suffix')):
    ilmn_suffix = x

nanopore_suffix = "nanopore.fq.gz"
if (x := validate_optional_key(config, 'NANOPORE_suffix')):
    nanopore_suffix = x

######################################
# Define some sample checking rules
######################################

def sampleid_valid_check(sampleid):
    if not sampleid.replace('_', '').isalnum():
        raise Exception(f"ERROR: {sampleid} contains non-alphanumeric or underscore characters.")
    else:
        return(True)

def sample_existence_check(top_raw_dir, sample_id, organisation_type = "folder", seq_platform = 'ILLUMINA'):
    if organisation_type == "folder":
        location = "{}/{}".format(top_raw_dir, sample_id)
        if os.path.exists(location):
            if seq_platform.upper() == 'ILLUMINA':
                if glob.glob("{}/{}/*[.-_]{}".format(top_raw_dir, sample_id, ilmn_suffix[0])):
                    return(True)
            elif seq_platform.upper() == 'NANOPORE':
                if glob.glob("{}/{}/*.{}".format(top_raw_dir, sample_id, nanopore_suffix)):
                    return(True)
    else:
        sample_pattern = "{}/{}[.-_]{}".format(top_raw_dir, sample_id, ilmn_suffix[0])
        if glob.glob(sample_pattern):
            return(True)
    return(False)

######################################
# Illumina-specific config parsing
######################################

# Make list of illumina samples, if ILLUMINA_SAMPLES in config
ilmn_samples = list()
ilmn_samples_organisation = "folder"
ilmn_runs_df = None
if (x := validate_optional_key(config, 'ILLUMINA_SAMPLES')):

    # Get relevant parameters
    raw_dir_ilmn         = validate_required_key(config, 'ILLUMINA_raw_reads_dir')
    perc_remaining_reads = validate_required_key(config, 'perc_remaining_reads')

    TRIMMOMATIC_threads  = validate_required_key(config, 'TRIMMOMATIC_threads')
    TRIMMOMATIC_memory   = validate_required_key(config, 'TRIMMOMATIC_memory')

    TRIMMOMATIC_adaptors = validate_required_key(config, 'TRIMMOMATIC_adaptors')
    if not os.path.exists(TRIMMOMATIC_adaptors):
        check_allowed_values('TRIMMOMATIC_adaptors', TRIMMOMATIC_adaptors, ('Skip', 'Quality'))

    if (TRIMMOMATIC_palindrome := validate_optional_key(config, 'TRIMMOMATIC_palindrome')):
        if not os.path.exists(TRIMMOMATIC_palindrome):
            raise Exception(f"ERROR in {config_path}: TRIMMOMATIC_palindrome file path does not exist.")

    if (TRIMMOMATIC_index_barcodes := validate_optional_key(config, 'TRIMMOMATIC_index_barcodes')):
        if TRIMMOMATIC_index_barcodes.lower() == 'infer':
            pass
        elif not os.path.exists(TRIMMOMATIC_index_barcodes):
            raise Exception(f"ERROR in {config_path}: TRIMMOMATIC_index_barcodes file path does not exist.")

    # Optional features with builtin defaults

    trimmomatic_simple_clip_threshold = 10
    if (y := validate_optional_key(config, 'TRIMMOMATIC_simple_clip_threshold')):
        trimmomatic_simple_clip_threshold = y

    multiplex_tech = 'TruSeq'
    if (y := validate_optional_key(config, 'ILLUMINA_MULTIPLEX_TECH')):
        multiplex_tech = y

    # column_name specified option
    if isinstance(x, str):
        ilmn_samples_organisation = "bulk"
        if os.path.exists(x) and os.path.isfile(x):
            # extra runs sheet
            print(f"MIntO uses {x} as sample list")
            col_name = "sample"
            ilmn_runs_df = pandas.read_table(x)
            for sampleid in ilmn_runs_df['sample'].unique():
                if sampleid_valid_check(sampleid):
                    for runid in ilmn_runs_df.loc[ilmn_runs_df['sample'] == sampleid]['run'].to_list():
                        #print(sampleid, runid)
                        if sample_existence_check(raw_dir_ilmn, runid, ilmn_samples_organisation):
                            if sampleid not in ilmn_samples:
                                ilmn_samples.append(sampleid)
                        else:
                            raise Exception(f"ERROR: {runid} not in bulk data folder {raw_dir_ilmn}")
        else:
            # column name in metadata sheet
            col_name = x
            md_df = pandas.read_table(metadata)
            if not col_name in md_df.columns:
                raise Exception(f"ERROR in {config_path}: column name specified for ILLUMINA_SAMPLES does not exist in metadata or runs sheet. Please fix!")
            sampleid_list = md_df[col_name].to_list()
            if sample_existence_check(raw_dir_ilmn, sampleid_list[0]):
                ilmn_samples_organisation = "folder"
            for sampleid in sampleid_list:
                if sampleid_valid_check(sampleid):
                    if sample_existence_check(raw_dir_ilmn, sampleid, ilmn_samples_organisation) and sampleid not in ilmn_samples:
                        ilmn_samples.append(sampleid)
                    else:
                        raise Exception(f"ERROR: {sampleid} not in raw data folder {raw_dir_ilmn} with suffix {ilmn_suffix[0]}")
    # listed samples
    else:
        for sampleid in x:
            if sampleid_valid_check(sampleid):
                if sample_existence_check(raw_dir_ilmn, sampleid):
                    ilmn_samples.append(sampleid)
                else:
                    raise Exception(f"ERROR: {sampleid} is not found under {raw_dir_ilmn} with suffix {ilmn_suffix[0]}.")

    # Register and process composite samples

    # Make list of merged illumina samples, if ILLUMINA_MERGE_SAMPLES in config
    merged_illumina_samples = dict()
    if (x := validate_optional_key(config, 'ILLUMINA_MERGE_SAMPLES')):
        for m in x:
            merged_illumina_samples.append(m)

    # If it's composite sample, then don't need to see them until it gets merged later
    x = ilmn_samples
    for i in x:
        if i in merged_illumina_samples.keys():
            ilmn_samples.remove(i)

    ilmn_remove_merged_samples = False
    if (x := validate_optional_key(config, 'ILLUMINA_MERGE_SAMPLES_REMOVE_CONTRIBUTORS')):
        ilmn_remove_merged_samples = x

######################################
# Nanopore-specific config parsing
######################################

# Make list of nanopore samples, if NANOPORE_SAMPLES in config
nano_samples = list()
merged_nanopore_samples = dict()
nano_remove_merged_samples = False
if (x := validate_optional_key(config, 'NANOPORE_SAMPLES')):
    raw_dir_nanopore = validate_required_key(config, 'NANOPORE_raw_reads_dir')
    for sampleid in x:
        if sampleid_valid_check(sampleid):
            if sample_existence_check(raw_dir_nanopore, sampleid, seq_platform='NANOPORE'):
                nano_samples.append(sampleid)
            else:
                raise Exception(f"ERROR: {sampleid} is not found under {raw_dir_nanopore} with suffix {nanopore_suffix}.")

    # Register and process composite samples

    # Make list of merged nanopore samples, if NANOPORE_MERGE_SAMPLES in config
    if (x := validate_optional_key(config, 'NANOPORE_MERGE_SAMPLES')):
        merged_nanopore_samples = x

    # If it's composite sample, then don't need to see them until it gets merged later
    x = nano_samples
    for i in x:
        if i in merged_nanopore_samples.keys():
            nano_samples.remove(i)

    if (x := validate_optional_key(config, 'NANOPORE_MERGE_SAMPLES_REMOVE_CONTRIBUTORS')):
        nano_remove_merged_samples = x

######################################
# Define all the outputs needed by target 'all'
######################################

def qc1_multiqc_output():
    result = []
    if len(ilmn_samples) > 0:
        result = expand("{wd}/output/1-0-qc/{omics}_multiqc.html",
                        wd = working_dir,
                        omics = omics)
    return(result)

def qc1_read_length_cutoff_output():
    result = []
    if len(ilmn_samples) > 0:
        result = expand("{wd}/output/1-trimmed/{omics}_cumulative_read_length_cutoff.pdf",
                        wd = working_dir,
                        omics = omics)
    return(result)

def qc1_config_yml_output():
    result = expand("{wd}/{omics}/QC_2.yaml",
                    wd = working_dir,
                    omics = omics)
    return(result)

def qc1_nanopore_output():
    result = []
    if len(nano_samples) > 0:
        result = expand("{wd}/{omics}/1-trimmed/{sample}/{sample}.nanopore.fq.gz",
                        wd = working_dir,
                        omics = omics,
                        sample = nano_samples)
    return(result)

rule all:
    input:
        qc1_nanopore_output(),
        qc1_multiqc_output(),
        qc1_read_length_cutoff_output(),
        qc1_config_yml_output(),
        print_versions.get_version_output(snakefile_name)
    default_target: True


###############################################################################################
###############################################################################################
##
## ILLUMINA-specific rules
##
###############################################################################################
###############################################################################################

if len(ilmn_samples) > 0:

    ##########
    # Filter reads by quality
    ##########

    # Get a sorted list of runs for a sample

    def get_runs_for_sample(sample):
        runs = [sample]
        if ilmn_samples_organisation == "folder":
            sample_dir = '{raw_dir_ilmn}/{sample}'.format(raw_dir_ilmn=raw_dir_ilmn, sample=sample)
            runs = [ f.name.split(ilmn_suffix[0])[0][:-1] for f in os.scandir(sample_dir) if f.is_file() and f.name.endswith(ilmn_suffix[0]) ]
        elif ilmn_runs_df is not None:
            runs = ilmn_runs_df.loc[ilmn_runs_df['sample'] == sample]['run'].to_list()
        return(sorted(runs))

    # Get files for a run

    def get_raw_reads_for_sample_run(wildcards):
        prefix = '{raw_dir_ilmn}/{run}'.format(raw_dir_ilmn=raw_dir_ilmn, run=wildcards.run)
        if ilmn_samples_organisation == "folder":
            prefix = '{raw_dir_ilmn}/{sample}/{run}'.format(raw_dir_ilmn=raw_dir_ilmn, sample=wildcards.sample, run=wildcards.run)
        raw_sample_run = {}
        for i, k in enumerate(['read_fw', 'read_rv']):
            raw_sample_run[k] = glob.glob("{}[.-_]{}".format(prefix, ilmn_suffix[i]))[0]
        return raw_sample_run

    # Get file for infering index

    def get_example_to_infer_index(sample):
        runs = get_runs_for_sample(sample)
        if ilmn_samples_organisation == "folder":
            sample_dir = '{raw_dir_ilmn}/{sample}'.format(raw_dir_ilmn=raw_dir_ilmn, sample=sample)
            example_fqs = [ os.path.normpath(f) for f in os.scandir(sample_dir) if f.is_file() and f.name.endswith(ilmn_suffix[0]) ]
        else:
            example_fqs = glob.glob("{}/{}[.-_]{}".format(raw_dir_ilmn, runs[0], ilmn_suffix[0]))
        return(example_fqs[0])

    ##########
    # FastQC of the initial raw data
    ##########

    rule initial_fastqc:
        input:
            unpack(get_raw_reads_for_sample_run)
        output:
            flag="{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastqc.done",
        params:
            outdir="{wd}/{omics}/1-0-qc/{sample}"
        log:
            "{wd}/logs/{omics}/1-0-qc/{sample}_{run}_qc1_fastqc.log"
        resources:
            mem=4
        threads:
            2
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            time (
                fastqc --noextract -f fastq -t {threads} -q -o {params.outdir} {input.read_fw} {input.read_rv} && \
            echo "Fastqc done" > {output.flag}
            ) >& {log}
            """

    rule qc1_fake_move_for_multiqc:
        input:
            fastqc_html=lambda wildcards: expand("{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastqc.done",
                                                wd=wildcards.wd,
                                                omics=wildcards.omics,
                                                sample=wildcards.sample,
                                                run=get_runs_for_sample(wildcards.sample)
                                                )
        output:
            temp("{wd}/{omics}/1-0-qc/{sample}_fake.moved")
        threads:
            2
        shell:
            """
            echo "{wildcards.sample} moved" > {output}
            """

    rule qc1_create_multiqc:
        input:
            fqc_flag=lambda wildcards: expand("{wd}/{omics}/1-0-qc/{sample}_fake.moved",
                                                wd=wildcards.wd,
                                                omics=wildcards.omics,
                                                sample=ilmn_samples
                                                )
        output:
            "{wd}/output/1-0-qc/{omics}_multiqc.html"
        params:
            indir="{wd}/{omics}/1-0-qc",
            outdir="{wd}/output/1-0-qc"
        log:
            "{wd}/logs/{omics}/1-0-qc/{omics}_multiqc.log"
        threads:
            2
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            time (
                multiqc --filename {wildcards.omics}_multiqc.html --outdir {params.outdir} -d --zip-data-dir -q --no-ansi --interactive {params.indir}
            ) >& {log}
            """

    ##########
    # Index barcodes should be used for adaptor trimming
    ##########

    if TRIMMOMATIC_index_barcodes:

        # Create a file with fwd/rev index sequences in the format
        #   GTTACGGA+CTTGGCTA
        # If there is a global index file for all samples given in config['TRIMMOMATIC_index_barcodes'], use it
        # If it does not exist, then infer from fastq file

        ruleorder: qc1_make_custom_adapter_file_with_palindrome > qc1_make_custom_adapter_file

        if TRIMMOMATIC_index_barcodes.lower() == "infer":
            # Infer from fastq file
            rule qc1_get_index_barcode:
                localrule: True
                input:
                    read_fw=lambda wildcards: get_example_to_infer_index(wildcards.sample),
                output:
                    barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
                shell:
                    """
                    outdir=$(dirname {output.barcodes})
                    mkdir -p $outdir
                    # Get the most frequent barcode-pair
                    # I am invoking a shell within, otherwise snakemake's strict mode catches SIG_PIPE from head to gzip and fails
                    barcodes=$(sh -c 'gzip -cd {input.read_fw} | grep "^@" | head -10000 | cut -f2 -d" " | cut -f4 -d":" | sort | uniq -c | sed "s/^ \+//" | tr -s " " " " | sort -k1,1nr -t" " | head -1 | cut -f2 -d" "')
                    # Write barcodes for this sample
                    echo $barcodes > {output.barcodes}
                    """

        else:
            # From global index file
            rule qc1_get_index_barcode:
                localrule: True
                input:
                    barcodes=TRIMMOMATIC_index_barcodes,
                output:
                    barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
                shell:
                    """
                    outdir=$(dirname {output.barcodes})
                    mkdir -p $outdir
                    # Get barcodes for this sample
                    barcodes=$(grep "^"{wildcards.sample}"[[:space:]]" {input.barcodes} | cut -f2)
                    # Write barcodes for this sample
                    echo $barcodes > {output.barcodes}
                    """

        # Create a custom adapter file to be used in Trimmomatic, given the custom index file above
        rule qc1_make_custom_adapter_file_with_palindrome:
            localrule: True
            input:
                barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
                palindrome=TRIMMOMATIC_palindrome,
                template=TRIMMOMATIC_adaptors
            output:
                adapter="{wd}/{omics}/1-trimmed/{sample}/adapters.fa"
            params:
                multiplex=multiplex_tech
            conda:
                minto_dir + "/envs/MIntO_base.yml" # seqtk
            shell:
                """
                # Get the 1st part of the adapter pair
                barcode1=$(cat {input.barcodes} | sed "s/;/+/" | cut -f1 -d'+')
                # Get the 2nd part of adapter pair
                barcode2=$(cat {input.barcodes} | sed "s/;/+/" | cut -f2 -d'+')

                # For MGI, revcomp barcode2
                if [ "{params.multiplex}" == "MGIEasy" ]; then
                    barcode2=$(echo $barcode2 | tr 'ATGC' 'TACG' | rev)
                fi

                # Make custom adapters for this sample using its index sequences
                #  1. Make palindromes
                cat {input.template} | seqkit grep --quiet -n -f {input.palindrome} -w 1000 -o - | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" | seqtk seq -r - | sed "s^>.*/^>PrefixPE-Ad/^" | tr -s '12' '21' > {output.adapter}
                #  2. Make 5-prime adapters
                cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" >> {output.adapter}
                #  3. Make 3-prime adapters
                cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" | seqtk seq -r - | tr -s '12' '21' | sed "s^/^_rc/^" >> {output.adapter}
                """

        # Create a custom adapter file to be used in Trimmomatic, given the custom index file above
        rule qc1_make_custom_adapter_file:
            localrule: True
            input:
                barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
                template=TRIMMOMATIC_adaptors
            output:
                adapter="{wd}/{omics}/1-trimmed/{sample}/adapters.fa"
            params:
                multiplex=multiplex_tech
            conda:
                minto_dir + "/envs/MIntO_base.yml" # seqtk
            shell:
                """
                # Get the 1st part of the adapter pair
                barcode1=$(cat {input.barcodes} | sed "s/;/+/" | cut -f1 -d'+')
                # Get the 2nd part of adapter pair
                barcode2=$(cat {input.barcodes} | sed "s/;/+/" | cut -f2 -d'+')

                # For MGI, revcomp barcode2
                if [ "{params.multiplex}" == "MGIEasy" ]; then
                    barcode2=$(echo $barcode2 | tr 'ATGC' 'TACG' | rev)
                fi

                # Make custom adapters for this sample using its index sequences
                #  2. Make 5-prime adapters
                cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" > {output.adapter}
                #  3. Make 3-prime adapters
                cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" | seqtk seq -r - | tr -s '12' '21' | sed "s^/^_rc/^" >> {output.adapter}
                """

    ##########
    # No index barcodes but adaptor trimming needed
    ##########

    elif TRIMMOMATIC_adaptors != 'Skip':

        # Hardlink standard file if there is no index_barcode file
        rule qc1_copy_fixed_adapter_file:
            localrule: True
            input:
                template=TRIMMOMATIC_adaptors
            output:
                adapter="{wd}/{omics}/1-trimmed/{sample}/adapters.fa"
            shell:
                """
                ln --force {input.template} {output.adapter}
                """

    ##########
    # Adaptor-trimming + Quality-trimming or Quality-trimming only?
    ##########

    # If there is adapter file, then priority is to use it
    ruleorder: qc1_trim_quality_and_adapter > qc1_trim_quality

    if TRIMMOMATIC_adaptors == 'Skip':

        # Fake a trim
        rule qc1_trim_quality:
            localrule: True
            input:
                unpack(get_raw_reads_for_sample_run)
            output:
                pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.fq.gz",
                pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.fq.gz",
            log:
                "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc1_trim_quality.log"
            shell:
                """
                time (
                    ln --force {input.read_fw} {output.pairead1}
                    ln --force {input.read_rv} {output.pairead2}
                    echo "Skipped trimming and linked raw files to trimmed files."
                ) >& {log}
                """

    else:

        # Trim for quality with Trimmomatic
        rule qc1_trim_quality:
            localrule: True
            input:
                unpack(get_raw_reads_for_sample_run)
            output:
                pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.fq.gz",
                singleread1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.single.fq.gz",
                pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.fq.gz",
                singleread2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.single.fq.gz",
                summary="{wd}/{omics}/1-trimmed/{sample}/{run}.trim.summary"
            shadow:
                "minimal"
            log:
                "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc1_trim_quality.log"
            resources:
                mem=TRIMMOMATIC_memory
            threads:
                TRIMMOMATIC_threads
            conda:
                minto_dir + "/envs/MIntO_base.yml"#trimmomatic
            shell:
                """
                remote_dir=$(dirname {output.pairead1})
                time (
                    trimmomatic PE -threads {threads} \
                        -summary {output.summary} \
                        -phred33 \
                        {input.read_fw} {input.read_rv} \
                        {wildcards.run}.1.fq.gz {wildcards.run}.1.single.fq.gz \
                        {wildcards.run}.2.fq.gz {wildcards.run}.2.single.fq.gz \
                        TRAILING:20 LEADING:5 SLIDINGWINDOW:4:20 \
                    && rsync -a * $remote_dir/
                ) >& {log}
                """

    # Trim for quality and adapter with Trimmomatic
    rule qc1_trim_quality_and_adapter:
        input:
            unpack(get_raw_reads_for_sample_run),
            adapter='{wd}/{omics}/1-trimmed/{sample}/adapters.fa'
        output:
            pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.fq.gz",
            singleread1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.single.fq.gz",
            pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.fq.gz",
            singleread2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.single.fq.gz",
            summary="{wd}/{omics}/1-trimmed/{sample}/{run}.trim.summary"
        shadow:
            "minimal"
        params:
            simple_clip_threshold=trimmomatic_simple_clip_threshold
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc1_trim_quality.log"
        resources:
            mem=TRIMMOMATIC_memory
        threads:
            TRIMMOMATIC_threads
        conda:
            minto_dir + "/envs/MIntO_base.yml"#trimmomatic
        shell:
            """
            remote_dir=$(dirname {output.pairead1})
            time (
                trimmomatic PE -threads {threads} \
                    -summary {output.summary} \
                    -phred33 \
                    {input.read_fw} {input.read_rv} \
                    {wildcards.run}.1.fq.gz {wildcards.run}.1.single.fq.gz \
                    {wildcards.run}.2.fq.gz {wildcards.run}.2.single.fq.gz \
                    TRAILING:20 LEADING:5 SLIDINGWINDOW:4:20 \
                    ILLUMINACLIP:{input.adapter}:2:30:{params.simple_clip_threshold}:1:TRUE \
                && rsync -a * $remote_dir/
            ) >& {log}
            """

    rule qc1_check_read_length:
        input:
            pairead=lambda wildcards: expand("{wd}/{omics}/1-trimmed/{sample}/{run}.{group}.fq.gz",
                                                wd=wildcards.wd,
                                                omics=wildcards.omics,
                                                sample=wildcards.sample,
                                                run=get_runs_for_sample(wildcards.sample),
                                                group=wildcards.group
                                                )
        output:
            length="{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.read_length.txt"
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}.{group}.check_read_length.log"
        resources:
            mem=4
        threads:
            2
        shell:
            """
            time (
                sh -c 'gzip -cd {input.pairead} | awk -v f="{wildcards.sample}_{wildcards.group}" "{{if(NR%4==2) print length(\$1),f}}" | sort -n | uniq -c > {output.length}'
            ) >& {log}
            """

    rule qc1_check_read_length_merge:
        localrule: True
        input:
            length=lambda wildcards: expand("{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.read_length.txt", wd=wildcards.wd, omics=wildcards.omics, sample=ilmn_samples, group=['1', '2'])
        output:
            readlen_dist="{wd}/{omics}/1-trimmed/samples_read_length.txt",
        log:
            "{wd}/logs/{omics}/1-trimmed/qc1_check_read_length_merge.log"
        threads: 1
        shell:
            """
            time (
                cat {input.length} > {output.readlen_dist}
            ) >& {log}
            """

    rule qc1_cumulative_read_len_plot:
        localrule: True
        input:
            readlen_dist=rules.qc1_check_read_length_merge.output.readlen_dist
        output:
            plot="{wd}/output/1-trimmed/{omics}_cumulative_read_length_cutoff.pdf",
            cutoff_file="{wd}/{omics}/1-trimmed/QC_1_min_len_read_cutoff.txt"
        params:
            cutoff=perc_remaining_reads,
        log:
            "{wd}/logs/{omics}/1-trimmed/plot_cumulative_read_len.log"
        resources:
            mem=TRIMMOMATIC_memory
        conda:
            minto_dir + "/envs/r_pkgs.yml"
        shell:
            """
            time (
                Rscript {script_dir}/QC_cumulative_read_length_plot.R --input {input.readlen_dist} --frac {params.cutoff} --out_plot {output.plot} --out_cutoff {output.cutoff_file}
            ) >& {log}
            """

    rule qc2_config_yml_file_illumina:
        localrule: True
        input:
            cutoff_file=rules.qc1_cumulative_read_len_plot.output.cutoff_file
        output:
            config_file=temp("{wd}/{omics}/QC_2.yaml.illumina")
        params:
            sample_list = ilmn_samples,
            remove_merged_samples = ilmn_remove_merged_samples,
            merge_samples_directive = '\n'.join([" {} : {}".format(i, merged_illumina_samples[i]) for i in merged_illumina_samples.keys()])
        log:
            "{wd}/logs/{omics}/qc2_config_yml_file.illumina.log"
        shell:
            """
            cat > {output.config_file} <<___EOF___

#########################
# Read length filtering
#########################

ILLUMINA_READ_minlen: $(cat {input.cutoff_file})

######################
# Optionally, do you want to merge replicates or make pseudo samples
# E.g:
# ILLUMINA_MERGE_SAMPLES:
#  sample1: rep1a+rep1b+rep1c
#  sample2: rep2a+rep2b+rep2c
#
# The above directive will make 2 new composite or pseudo samples at the end of QC_2.
# Imagine you had triplicates for sample1 named as rep1a, rep1b and rep1c.
# And likewise for sample2. The directive above will:
#     - combine 3 samples (rep1a, rep1b and rep1c)into a new sample called 'sample1'.
#     - combine 3 samples (rep2a, rep2b and rep2c)into a new sample called 'sample2'.
# For all subsequent steps, namely:
#     - profiling (done by QC_2.smk),
#     - assembly (done by assembly.smk),
#     - binning (done by binning_preparation.smk and mags_generation.smk),
# you can just use 'sample2' instead of the replicates rep2a, rep2b and rep2c in the yaml files.
# Please note that METADATA file must have an entry for 'sample2' as well,
# otherwise QC_2 step will fail.
# Having extra entries in METADATA file does not affect you in any way.
# Therefore, it is safe to have metadata recorded for
# rep2a, rep2b, rep2c, sample2 from the beginning.
######################

ILLUMINA_MERGE_SAMPLES_REMOVE_CONTRIBUTORS: {params.remove_merged_samples}
ILLUMINA_MERGE_SAMPLES:
{params.merge_samples_directive}

# ILLUMINA_SAMPLES section:
# -----------------
# List of illumina samples that will be filtered by read length.
#
# E.g.:
# - S1
# - S2

ILLUMINA_SAMPLES:
$(for i in {params.sample_list}; do echo "- '$i'"; done)

___EOF___
            echo "ILLUMINA_SAMPLES: {params.sample_list}" >& {log}
            """

###############################################################################################
###############################################################################################
##
## NANOPORE-specific rules
##
###############################################################################################
###############################################################################################

if len(nano_samples) > 0:

    #################
    # Get HQ reads covering 95% bases
    #################

    rule filter_nanopore_reads:
        input:
            ont="{wd}/{omics}/0-raw/{sample}/{sample}.nanopore.fq.gz"
        output:
            summary="{wd}/{omics}/1-trimmed/{sample}/{sample}.NANOPORE.trim.summary",
            ont="{wd}/{omics}/1-trimmed/{sample}/{sample}.nanopore.fq.gz"
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_qc1_nanopore_trimlog.log"
        shadow:
            "minimal"
        threads: 4
        resources:
            mem=5
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            time (
                fastplong --disable_adapter_trimming --cut_front --cut_front_window_size 1 --cut_front_mean_quality 15 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 15 --thread {threads} -z 6 -i {input.ont} --stdout \
                        | tee >(seqkit stats --tabular --all -i {wildcards.sample} - > {output.summary}) \
                        | gzip > filtered.fq.gz
                rsync -a filtered.fq.gz {output.ont}
            ) >& {log}
            """


    rule qc2_config_yml_file_nanopore:
        localrule: True
        output:
            config_file=temp("{wd}/{omics}/QC_2.yaml.nanopore")
        params:
            sample_list = nano_samples,
            remove_merged_samples = nano_remove_merged_samples,
            merge_samples_directive = '\n'.join([" {} : {}".format(i, merged_nanopore_samples[i]) for i in merged_nanopore_samples.keys()])
        log:
            "{wd}/logs/{omics}/qc2_config_yml_file.nanopore.log"
        shell:
            """
            cat > {output.config_file} <<___EOF___
#########################
# Read length filtering
#########################

NANOPORE_READ_minlen: 500

######################
# Optionally, if you want to merge replicates or make pseudo samples
# E.g:
# NANOPORE_MERGE_SAMPLES:
#  sample1: rep1a+rep1b+rep1c
#  sample2: rep2a+rep2b+rep2c
#
# The above directive will make 2 new composite or pseudo samples at the end of QC_2.
# Imagine you had triplicates for sample1 named as rep1a, rep1b and rep1c.
# And likewise for sample2. The directive above will:
#     - combine 3 samples (rep1a, rep1b and rep1c)into a new sample called 'sample1'.
#     - combine 3 samples (rep2a, rep2b and rep2c)into a new sample called 'sample2'.
# For all subsequent steps, namely:
#     - profiling (done by QC_2.smk),
#     - assembly (done by assembly.smk),
#     - binning (done by binning_preparation.smk and mags_generation.smk),
# you can just use 'sample2' instead of the replicates rep2a, rep2b and rep2c in the yaml files.
# Please note that METADATA file must have an entry for 'sample2' as well,
# otherwise QC_2 step will fail.
# Having extra entries in METADATA file does not affect you in any way.
# Therefore, it is safe to have metadata recorded for
# rep2a, rep2b, rep2c, sample2 from the beginning.
######################

NANOPORE_MERGE_SAMPLES_REMOVE_CONTRIBUTORS: {params.remove_merged_samples}
NANOPORE_MERGE_SAMPLES:
{params.merge_samples_directive}

# NANOPORE_SAMPLES section:
# -----------------
# List of nanopore samples that will be filtered by read length.
#
# E.g.:
# - S1
# - S2

NANOPORE_SAMPLES:
$(for i in {params.sample_list}; do echo "- '$i'"; done)
___EOF___

            echo "NANOPORE_SAMPLES: {params.sample_list}" >& {log}
            """

###############################################################################################
###############################################################################################
##
## Generate core part of the yaml file
##
###############################################################################################
###############################################################################################

##########################################################################################################
# Generate configuration yml file for the next step in pre-processing of metaG and metaT data
##########################################################################################################

rule qc2_config_yml_core:
    localrule: True
    output:
        config_file=temp("{wd}/{omics}/QC_2.yaml.core")
    log:
        "{wd}/logs/{omics}/qc2_config_yml_file.core.log"
    shell:
        """
        # Generate part1 of config file
        cat > {output.config_file} <<___EOF___
######################
# General settings
######################

PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
minto_dir: {minto_dir}
METADATA: {metadata}

########################################
# Program settings
########################################

#########################
# Host genome filtering
#########################

# Path of the host genome file: <PATH_host_genome>/<NAME_host_genome>
# Inferred location for bwa-mem2/strobealign index files: <PATH_host_genome>/{{BWA|STROBEALIGN}}_index/<NAME_host_genome>.*
# If index files already exist, then they will be used directly.
# If not, a fasta file should exist exactly as: <PATH_host_genome>/<NAME_host_genome>
#   which will be used to build the index files in the inferred location above.

PATH_host_genome: None
NAME_host_genome: None

# Should host genome index files be cached in a local mirror?
# If yes, set the path here.
# If no,  set it to None

LOCAL_DATABASE_CACHE_DIR: None

# Which aligner or mapper to use for ILLUMINA: 'bwa' or 'strobealign' is supported

ILLUMINA_ALIGNER_type: bwa
ILLUMINA_ALIGNER_threads: 8

# Which aligner or mapper to use for NANOPORE: 'minimap2' is supported

NANOPORE_ALIGNER_type: minimap2
NANOPORE_ALIGNER_threads: 8
___EOF___

            # Generate part2 of config file
            if [ "{omics}" == "metaT" ]; then
                cat >> {output.config_file} <<___EOF___

#########################
# Ribosomal RNA depletion
#########################

sortmeRNA_threads: 8
sortmeRNA_memory: 10
sortmeRNA_db: {minto_dir}/data/rRNA_databases
sortmeRNA_db_idx: {minto_dir}/data/rRNA_databases/idx
___EOF___

            fi

            # Generate part3 of config file
            cat >> {output.config_file} <<___EOF___

##################################
# Assembly-free taxonomy profiling
##################################

# Following values for 'TAXA_profiler' are supported:
#    1. metaphlan - relative abundance using MetaPhlAn
#    2. motus_raw - read counts using mOTUs
#    3. motus_rel - relative abundance using mOTUs
# Comma-delimited combination of multiple options also supported
# Eg:
#    TAXA_profiler: metaphlan,motus_rel
TAXA_threads: 8
TAXA_memory: 15
TAXA_profiler: motus_rel,metaphlan
metaphlan_version: 4.2.2
motus_version: 3.1.0

#########################
# K-mer based comparison
#########################

# FracMinHash comparisons by sourmash
# SOURMASH_min_abund - Minimum count of each k-mer for filtering the sketch (integer)
# SOURMASH_max_abund - Maximum count of each k-mer for filtering the sketch (integer)
# SOURMASH_cutoff    - Dissimilarity cutoff for subclusters via hierarchical clustering

SOURMASH_min_abund: 2
SOURMASH_max_abund: 1000
SOURMASH_cutoff: 0.40

#####################
# Analysis parameters
#####################

# MAIN_factor  - the main factor in the metadata file to differentiate in visualization (using color)
# PLOT_factor2 - the second factor in the metadata file to differentiate in visualization (using shape)
# PLOT_time    - name of the factor in the metadata file denoting time (e.g. hour, day)

MAIN_factor:
PLOT_factor2:
PLOT_time:

#####################
# Co-assembly grouping
#####################

# COAS_factor - The factor/attribute to use to group samples for co-assembly

COAS_factor:

######################
# Input data
######################

___EOF___
        """

###############################################################################################
###############################################################################################
##
## Combine tech-specific yaml files
##
###############################################################################################
###############################################################################################

def get_qc2_config_pieces(wildcards):
    wd = wildcards.wd
    omics = wildcards.omics

    results = list()
    results.append(f"{wd}/{omics}/QC_2.yaml.core")
    if len(ilmn_samples) > 0:
        results.append(f"{wd}/{omics}/QC_2.yaml.illumina")
    if len(nano_samples) > 0:
        results.append(f"{wd}/{omics}/QC_2.yaml.nanopore")
    return(results)

rule qc2_filter_config_yml_assembly:
    localrule: True
    input:
        pieces = get_qc2_config_pieces,
    output:
        config_file="{wd}/{omics}/QC_2.yaml"
    resources:
        mem=2
    threads: 1
    log:
        "{wd}/logs/{omics}/config_yml_QC_2.log"
    shell:
        """
        cat {input.pieces} > {output}
        """
