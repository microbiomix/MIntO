#!/usr/bin/env python

'''
Assembling metagenomes from combinations of illumina/bgi-seq and nanopore sequencing data.

Authors: Carmen Saenz, Mani Arumugam
'''

import os

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/locations.smk'
include: 'include/fasta_bam_helpers.smk'
include: 'include/resources.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule assembly_base from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# Variables from configuration yaml file

ilmn_samples            = list()
nanopore_samples        = list()
merged_illumina_samples = dict()

# Common variables

spades_script = 'spades.py' # from conda environment

##############################################
# Register composite samples
##############################################

# Make list of illumina coassemblies, if ILLUMINA_MERGE_SAMPLES in config
if (x := validate_optional_key(config, 'ILLUMINA_MERGE_SAMPLES')):
    for m in x:
        merged_illumina_samples.append(m)

##############################################
# Get sample list
##############################################

# Make list of illumina samples, if ILLUMINA_SAMPLES in config
if (x := validate_optional_key(config, 'ILLUMINA_SAMPLES')):
    check_input_directory(x, locations = ['6-corrected', '5-corrected-runs', get_qc2_output_location(omics)])
    ilmn_samples = x

    # If it's composite sample, then don't need to see them until it gets merged later
    for i in x:
        if i in merged_illumina_samples.keys():
            ilmn_samples.remove(i)


# Make list of nanopore samples, if NANOPORE_SAMPLES in config
if (x := validate_optional_key(config, 'NANOPORE_SAMPLES')):
    check_input_directory(x, locations = ['6-corrected', get_qc2_output_location(omics)])
    nanopore_samples = x

###############################
# Hybrid-/co-assemblies
###############################

# Make list of nanopore-illumina hybrid assemblies, if HYBRID in config
hybrid_assemblies = list()
if (x := validate_optional_key(config, 'HYBRID')):
    print(f"Found HYBRID in {config_path}. Enabling hybrid assembly.")
    for nano in x:
        for ilmn in x[nano].split("+"):
            hybrid_assemblies.append(nano+"-"+ilmn)

# Make list of illumina coassemblies, if enable_COASSEMBLY is set to "yes" in config
co_assemblies = dict()
if validate_optional_key(config, 'enable_COASSEMBLY'):
    if (x := validate_optional_key(config, 'COASSEMBLY')):
        print(f"Found COASSEMBLY in {config_path}. Enabling co-assembly.")
        co_assemblies = x

###############################
# ILLUMINA-specific parameters
###############################

if len(ilmn_samples) > 0:

    ###############################
    # MetaSPAdes parameters
    ###############################

    METASPADES_qoffset = validate_required_key(config, 'METASPADES_qoffset')
    check_allowed_values('METASPADES_qoffset', METASPADES_qoffset, ('auto', '33', '64'))

    METASPADES_threads = validate_required_key(config, 'METASPADES_threads')

    METASPADES_illumina_max_k = validate_required_key(config, 'METASPADES_illumina_max_k')
    check_number_is_odd('METASPADES_illumina_max_k', METASPADES_illumina_max_k)

    METASPADES_hybrid_max_k = validate_optional_key(config, 'METASPADES_hybrid_max_k')
    if METASPADES_hybrid_max_k is not None:
        check_number_is_odd('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k)

    # Figure out SPAdes build
    # and verify kmer max

    if (x := validate_optional_key(config, 'METASPADES_custom_build')):
        spades_script = x # custom built, e.g. for higher K
        check_number_is_between('METASPADES_illumina_max_k', METASPADES_illumina_max_k, 19, 300)
        if METASPADES_hybrid_max_k is not None:
            check_number_is_between('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k, 19, 300)
    else:
        check_number_is_between('METASPADES_illumina_max_k', METASPADES_illumina_max_k, 19, 128)
        if METASPADES_hybrid_max_k is not None:
            check_number_is_between('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k, 19, 128)

    MEGAHIT_threads = validate_required_key(config, 'MEGAHIT_threads')
    MEGAHIT_memory  = validate_required_key(config, 'MEGAHIT_memory')

    ###############################
    # MEGAHIT parameter sets
    ###############################

    MEGAHIT_presets = list()
    mega_k_list = []

    # Check for MEGAHIT_presets
    if (x := validate_optional_key(config, 'MEGAHIT_presets')):
        MEGAHIT_presets = x

    # Also check for MEGAHIT_custom
    if (x := validate_optional_key(config, 'MEGAHIT_custom')):
        # if the custom k-s are set, that should be added to the MEGAHIT assembly types
        if isinstance(x, str):
            x = [x]
        for i, k_list in enumerate(x):
            if k_list and not k_list.isspace():
                MEGAHIT_presets.append(f'meta-custom-{i+1}')
                mega_k_list.append(k_list)

    if co_assemblies and not MEGAHIT_presets:
        raise Exception(f"ERROR in {config_path}: MEGAHIT_presets list of MEGAHIT parameters to run per co-assembly is empty")

###############################
# NANOPORE-specific parameters
###############################

if len(nanopore_samples) > 0:

    ###############################
    # MetaFlye parameter sets
    ###############################

    # Check for METAFLYE_presets

    METAFLYE_presets = validate_required_key(config, 'METAFLYE_presets')
    if (x := validate_optional_key(config, 'MEDAKA_INFERENCE_MODEL')):
        MEDAKA_INFERENCE_MODEL = x

# Define all the outputs needed by target 'all'

def illumina_single_assembly_output():
    if len(ilmn_samples) == 0:
        return(list())

    result = expand("{wd}/{omics}/7-assembly/{sample}/{kmer_dir}/{sample}.{sequence}.fasta",
                    wd = working_dir,
                    omics = omics,
                    sample = ilmn_samples,
                    kmer_dir = "k21-" + str(METASPADES_illumina_max_k),
                    sequence = ["contigs", "scaffolds"])
    return(result)

def illumina_co_assembly_output():
    if len(co_assemblies) == 0:
        return(list())

    result = expand("{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta",
                    wd = working_dir,
                    omics = omics,
                    coassembly = co_assemblies,
                    assembly_preset = MEGAHIT_presets)
    return(result)

def nanopore_single_assembly_output():
    if len(nanopore_samples) == 0:
        return(list())

    result = expand("{wd}/{omics}/7-assembly/{sample}/{assembly_preset}/{sample}.assembly.{seq_type}.fasta",
                    wd = working_dir,
                    omics = omics,
                    sample = nanopore_samples,
                    seq_type = ['nodes', 'edges'],
                    assembly_preset = METAFLYE_presets)
    return(result)

def hybrid_assembly_output():
    if len(hybrid_assemblies) == 0:
        return(list())

    result = expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.{sequence}.fasta",
                    wd = working_dir,
                    omics = omics,
                    assembly = hybrid_assemblies,
                    kmer_dir = "k21-" + str(METASPADES_hybrid_max_k),
                    sequence = ["contigs", "scaffolds"])
    return(result)

# Figure out kmer details

def get_metaspades_kmer_option(kk):
    kmers = [21] # start at 21
    kmers.extend(list(range(33, kk, 22))) # add more kmers with 22 increment while < max_k
    kmers.extend([kk]) # add max_k
    kmer_option = ','.join([str(k) for k in kmers])
    return(kmer_option)

rule all:
    input:
        illumina_single_assembly_output(),
        illumina_co_assembly_output(),
        nanopore_single_assembly_output(),
        hybrid_assembly_output(),
        print_versions.get_version_output(snakefile_name)
    default_target: True

###############################
# Common functions/rules
###############################

def get_corrected_runs_for_sample(wildcards):
    if wildcards.pair == 'nanopore':
        files = expand("{wd}/{omics}/4-hostfree/{nanopore}/{run}.{pair}.fq.gz",
                                    wd = wildcards.wd,
                                    omics = wildcards.omics,
                                    nanopore = wildcards.nanopore,
                                    run = get_runs_for_sample(wildcards.wd, wildcards.omics, wildcards.nanopore, caller='assembly', seq_platform='NANOPORE'),
                                    pair = wildcards.pair)
    else:
        if wildcards.illumina in merged_illumina_samples:
            files = list()
            for r in merged_illumina_samples[wildcards.illumina].split('+'):
                rep = r.strip()
                files += expand("{wd}/{omics}/5-corrected-runs/{rep}/{run}.{pair}.fq.gz",
                                        wd = wildcards.wd,
                                        omics = wildcards.omics,
                                        rep = rep,
                                        run = get_runs_for_sample(wildcards.wd, wildcards.omics, rep, caller='assembly'),
                                        pair = wildcards.pair)
        else:
            files = expand("{wd}/{omics}/5-corrected-runs/{illumina}/{run}.{pair}.fq.gz",
                                        wd = wildcards.wd,
                                        omics = wildcards.omics,
                                        illumina = wildcards.illumina,
                                        run = get_runs_for_sample(wildcards.wd, wildcards.omics, wildcards.illumina, caller='assembly', seq_platform='ILLUMINA'),
                                        pair = wildcards.pair)
    return(files)

rule merge_runs_base:
    localrule: True
    input:
        files=get_corrected_runs_for_sample
    output:
        combined="{wd}/{omics}/6-corrected/{illumina}/{illumina}.{pair}.fq.gz.impossible"
    wildcard_constraints:
        pair="1|2|nanopore"
    shadow:
        "minimal"
    params:
        num_runs = lambda wildcards, input: len(input.files)
    shell:
        """
        if (( {params.num_runs} > 1 )); then
            cat {input} > combined.fq.gz
            rsync -a combined.fq.gz {output}
        else
            cd $(dirname {output})
            ln --force {input[0]} {output}
        fi
        """

###############################
# NANOPORE-specific rules
###############################

if len(ilmn_samples) > 0:

    ###############################################################################################
    # Correct Illumina reads using SPAdes' spadeshammer
    ###############################################################################################

    def get_hq_fastq_files_illumina(wildcards):
        return(expand("{wd}/{omics}/{location}/{illumina}/{run}.{pair}.fq.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    location=get_qc2_output_location(wildcards.omics),
                    illumina = wildcards.illumina,
                    run = wildcards.run,
                    pair = [1, 2]))

    rule correct_spadeshammer:
        input:
            reads=get_hq_fastq_files_illumina
        output:
            fwd=temp("{wd}/{omics}/5-corrected-runs/{illumina}/{run}.1.fq.gz"),
            rev=temp("{wd}/{omics}/5-corrected-runs/{illumina}/{run}.2.fq.gz"),
        shadow:
            "minimal"
        params:
            qoffset = METASPADES_qoffset
        resources:
            mem = lambda wildcards, input, attempt: int(26 + 52*get_file_size_gb(input.reads[0]) + 50*(attempt-1))
        log:
            "{wd}/logs/{omics}/5-corrected-runs/{illumina}/{run}_spadeshammer.log"
        threads:
            METASPADES_threads
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            mkdir -p $(dirname {output.fwd})
            time (
                rsync --copy-links {input.reads[0]} {wildcards.run}.1.fq.gz
                rsync --copy-links {input.reads[1]} {wildcards.run}.2.fq.gz
                {spades_script} --only-error-correction -1 {wildcards.run}.1.fq.gz -2 {wildcards.run}.2.fq.gz -t {threads} -m {resources.mem} -o {wildcards.run} --phred-offset {params.qoffset}
                rsync -a {wildcards.run}/corrected/{wildcards.run}.1.fq00.0_0.cor.fastq.gz {output.fwd}; rsync -a {wildcards.run}/corrected/{wildcards.run}.2.fq00.0_0.cor.fastq.gz {output.rev}
            ) >& {log}
            """

    use rule merge_runs_base as merge_runs_illumina with:
        output:
            combined="{wd}/{omics}/6-corrected/{illumina}/{illumina}.{pair}.fq.gz"
        wildcard_constraints:
            pair="1|2"

    ruleorder: hybrid_assembly_metaspades > illumina_assembly_metaspades

    ###############################################################################################
    ########  Individual assembly of illumina samples
    ###############################################################################################
    rule illumina_assembly_metaspades:
        input:
            fwd="{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz",
            rev="{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz",
        output:
            cont_fa    = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/contigs.fasta",
            cont_pth   = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/contigs.paths",
            scaf_fa    = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/scaffolds.fasta",
            scaf_pth   = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/scaffolds.paths",
            scaf_gfa   = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/assembly_graph_with_scaffolds.gfa.gz",
            asm_log    = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/spades.log",
            asm_params = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/params.txt",
        shadow:
            "minimal"
        params:
            qoffset        = METASPADES_qoffset,
            asm_mode       = "--meta",
            kmer_option    = lambda wildcards: get_metaspades_kmer_option(int(wildcards.maxk)),
            nanopore_input = ""
        resources:
            mem = lambda wildcards, input, attempt: int(26 + 8.8*get_file_size_gb(input.fwd) + 20*(attempt-1))
        log:
            "{wd}/logs/{omics}/7-assembly/{illumina}/k21-{maxk}/{illumina}_metaspades.log"
        threads:
            METASPADES_threads
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            time (
                rsync --copy-links {input.fwd} {wildcards.illumina}.1.fq.gz
                rsync --copy-links {input.rev} {wildcards.illumina}.2.fq.gz
                {spades_script} {params.asm_mode} --only-assembler -1 {wildcards.illumina}.1.fq.gz -2 {wildcards.illumina}.2.fq.gz {params.nanopore_input} -t {threads} -m {resources.mem} -o outdir --tmp-dir tmp --phred-offset {params.qoffset} -k {params.kmer_option}
                rsync -a outdir/contigs.fasta {output.cont_fa}
                rsync -a outdir/contigs.paths {output.cont_pth}
                rsync -a outdir/scaffolds.fasta {output.scaf_fa}
                rsync -a outdir/scaffolds.paths {output.scaf_pth}
                rsync -a outdir/spades.log {output.asm_log}
                rsync -a outdir/params.txt {output.asm_params}
                gzip -c outdir/assembly_graph_with_scaffolds.gfa > {output.scaf_gfa}
            ) >& {log}
            """

    ###############################################################################################
    # Hybrid assembly of illumina and nanopore samples
    #  - inherits the normal illumina assembly rule above with --nanopore argument
    ###############################################################################################
    use rule illumina_assembly_metaspades as hybrid_assembly_metaspades with:
        input:
            fwd="{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz",
            rev="{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz",
            ont="{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz"
        output:
            cont_fa    = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/contigs.fasta",
            cont_pth   = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/contigs.paths",
            scaf_fa    = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/scaffolds.fasta",
            scaf_pth   = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/scaffolds.paths",
            scaf_gfa   = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/assembly_graph_with_scaffolds.gfa.gz",
            asm_log    = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/spades.log",
            asm_params = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/params.txt",
        params:
            qoffset        = METASPADES_qoffset,
            asm_mode       = "--meta",
            kmer_option    = lambda wildcards: get_metaspades_kmer_option(int(wildcards.maxk)),
            nanopore_input = lambda wildcards, input: "--nanopore {}".format(input.ont)
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/{nanopore}-{illumina}_metaspades.log"

    ###############################################################################################
    # This checks whether the first and last k characters in the contig are the same.
    # If yes, it removes the terminal k characters and appends '_circularA' to the fasta header.
    # It also updates 'length' in the fasta header.
    ###############################################################################################
    rule mark_circular_metaspades_contigs:
        localrule: True
        input:
            "{wd}/{omics}/7-assembly/{sample}/k21-{maxk}/{which}.fasta"
        output:
            "{wd}/{omics}/7-assembly/{sample}/k21-{maxk}/{sample}.{which}.fasta"
        params:
            kmer = lambda wildcards: int(wildcards.maxk)
        run:
            import re
            regex = re.compile('length_[0-9]+_')
            # Open the output file
            with open(output[0], 'w') as out:
                # Go through the fasta file
                fiter = fasta_iter(input[0])
                for entry in fiter:
                    header, seq = entry
                    if len(seq) > params.kmer:
                        if seq[0:params.kmer] == seq[-params.kmer:]:
                            seq = seq[0:-params.kmer]
                            header = regex.sub("length_%s_" % len(seq), header) + '_circularA'
                    out.write(">MetaSPAdes.k21-{maxk}.{sample}_{header}\n".format(maxk=wildcards.maxk, sample=wildcards.sample, header=header))
                    out.write(seq+"\n")

    ###############################################################################################
    ########  Co-assembly
    # This starts with 11G per sample, but if that fails, it increases by 5G per sample per repeated attempt
    ###############################################################################################

    def get_megahit_parameters(wildcards, kk, k_list):
        if (wildcards.assembly_preset == 'rnaspades'):
            # We assume kk is about 1/2 max-length - see rnaSPAdes recommendation
            # We will then make kmers=[0.33max, ..., 0.5max] with step 22
            min_k = 2*int(kk/3)+1 # start at odd number near 1/3 max-length
            kmers = list(range(min_k, kk-1, 22)) # add more kmers with 22 increment while < max_k
            kmers.extend([kk])
            kmer_option = ','.join([str(k) for k in kmers])
            return "--k-list {}".format(kmer_option)
        elif (wildcards.assembly_preset.startswith('meta-custom')):
            k_list_i = int(wildcards.assembly_preset.rsplit("-", 1)[-1])
            return "--k-list {}".format(k_list[k_list_i - 1])
        elif (wildcards.assembly_preset in ['meta-large', 'meta-sensitive']):
            return "--presets {}".format(wildcards.assembly_preset)
        else:
            return ""

    rule coassembly_megahit:
        input:
            fwd=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz', wd=wildcards.wd, omics=wildcards.omics, illumina=co_assemblies[wildcards.coassembly].split('+')),
            rev=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz', wd=wildcards.wd, omics=wildcards.omics, illumina=co_assemblies[wildcards.coassembly].split('+'))
        output:
            cont_fa    = "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/final.contigs.fa",
            asm_log    = "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/log.txt",
            asm_params = "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/options.json",
        shadow:
            "minimal"
        params:
            fwd_reads=lambda wildcards, input: ",".join(input.fwd),
            rev_reads=lambda wildcards, input: ",".join(input.rev),
            asm_params=lambda wildcards: get_megahit_parameters(wildcards, METASPADES_illumina_max_k, mega_k_list),
        resources:
            mem       = lambda wildcards, input, attempt: min(900, len(input.fwd)*(MEGAHIT_memory+6*(attempt-1))),
            mem_bytes = lambda wildcards, input, attempt: min(900, len(input.fwd)*(MEGAHIT_memory+6*(attempt-1)))*1024*1024*1024
        log:
            "{wd}/logs/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}_{assembly_preset}_coassembly_megahit.log"
        threads:
            MEGAHIT_threads
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            # Don't create the --out-dir directory as MEGAHIT wants it to not exist before
            time (
                megahit -1 {params.fwd_reads} -2 {params.rev_reads} -t {threads} -m {resources.mem_bytes} --out-dir outdir {params.asm_params}
            ) >& {log}
            rsync -a outdir/final.contigs.fa {output.cont_fa}
            rsync -a outdir/log {output.asm_log}
            rsync -a outdir/options.json {output.asm_params}
            """

    # Rename MEGAHIT contigs
    rule rename_megahit_contigs:
        localrule: True
        input:
            "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/final.contigs.fa"
        output:
            "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta"
        wildcard_constraints:
            assembly_preset = '|'.join(MEGAHIT_presets)
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            r"""
            perl -ne 's/^>k(\d+)_(\d+) (.*)len=(\d+)/>MEGAHIT.{wildcards.assembly_preset}.{wildcards.coassembly}_NODE_$2_length_$4_k_$1/ if m/^>/; print $_;' < {input} > {output}
            """

###############################
# NANOPORE-specific rules
###############################

if len(nanopore_samples) > 0:

    use rule merge_runs_base as merge_runs_nanopore with:
        output:
            combined="{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.{pair}.fq.gz"
        wildcard_constraints:
            pair="nanopore"

    ###############################################################################################
    # Assemble nanopore sequences using flye
    ###############################################################################################
    rule nanopore_assembly_metaflye:
        input:
            ont   = "{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz"
        output:
            fasta = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly.fasta",
            gfa   = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_graph.gfa",
            gv    = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_graph.gv",
            info  = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_info.txt"
        shadow:
            "minimal"
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}_{assembly_preset}_metaflye.log"
        resources:
            mem = lambda wildcards, attempt: 30*attempt
        threads: 16
        params:
            options = lambda wildcards: METAFLYE_presets[wildcards.assembly_preset]
        conda:
            minto_dir + "/envs/MIntO_base.yml"
        shell:
            """
            mkdir -p $(dirname {output[0]})
            time (
                flye --nano-hq {input} --out-dir out --threads {threads} --meta {params.options}
                rsync -a out/assembly.fasta {output.fasta}
                rsync -a out/assembly_info.txt {output.info}
                rsync -a out/assembly_graph.gfa {output.gfa}
                rsync -a out/assembly_graph.gv {output.gv}
            ) >& {log}
            """

    ###############################################################################################
    # Dnaapler reorientation for mixed metagenomic FASTA
    # "dnaapler all" sweeps for chromosomes, plasmids, and phages simultaneously.
    # We ignore sequences not flagged as circular by Flye
    ###############################################################################################
    rule dnaapler_reorient_nodes:
        input:
            contigs = rules.nanopore_assembly_metaflye.output.fasta,
            info    = rules.nanopore_assembly_metaflye.output.info
        output:
            fasta   = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly.dnaapler.fasta"
        shadow:
            "minimal"
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/dnaapler_fasta.log"
        threads: 8
        conda:
            minto_dir + "/envs/ont_polishing.yml"
        shell:
            """

            awk '$4 == "N"' {input.info} | cut -f 1 > ignore.list

            time (
                dnaapler all -i {input.contigs} -o out --ignore ignore.list -t {threads} -f
            ) >& {log}

            # Extract Dnaapler's aggregated output file, fallback safely if none were rotated
            if [ -f "out/dnaapler_reoriented.fasta" ]; then
                rsync -a "out/dnaapler_reoriented.fasta" {output.fasta}
            else
                ln --force {input.contigs} {output.fasta}
            fi
            """

    ###############################################################################################
    # Dnaapler reorientation for mixed metagenomic gfa
    # "dnaapler all" sweeps for chromosomes, plasmids, and phages simultaneously.
    # Since gfa mode only processes circular edges, we don't need an ignore list here.
    ###############################################################################################
    rule dnaapler_reorient_edges:
        input:
            gfa = rules.nanopore_assembly_metaflye.output.gfa,
        output:
            gfa = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_graph.dnaapler.gfa"
        shadow:
            "minimal"
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/dnaapler_gfa.log"
        threads: 8
        conda:
            minto_dir + "/envs/ont_polishing.yml"
        shell:
            """
            set -euo pipefail
            (
                set +e
                time dnaapler all -i {input.gfa} -o out -t {threads} -f
                dnaapler_status=$?
                set -e

                echo "dnaapler exit status: $dnaapler_status"

                # Extract Dnaapler's aggregated output file, fallback safely if none were rotated
                if [ -f "out/dnaapler_reoriented.gfa" ]; then
                    echo "dnaapler produced reoriented GFA; using it."
                    rsync -a "out/dnaapler_reoriented.gfa" {output.gfa}
                else
                    echo "dnaapler did not produce reoriented GFA; using original GFA."
                    ln --force {input.gfa} {output.gfa}
                fi
            ) >& {log}
            """

    ###############################################################################################
    # Create an edge-fasta from the dnaapler reoriented gfa file
    # For reoriented edges, it adds "rotated=True rotated_gene=dnaA" dnaaper-style fasta description
    # Also write a assembly_graph_info.txt where the first four columns are just like
    # Flye's assembly_info.txt on nodes
    ###############################################################################################

    rule extract_gfa_edge_info_fasta:
        input:
            gfa   = rules.dnaapler_reorient_edges.output.gfa,
        output:
            fasta = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_graph.dnaapler.fasta",
            info  = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_graph_info.txt"
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/extract_gfa_edge_info_fasta.log"
        resources:
            mem = 5
        run:
            import textwrap

            edge_order = []
            seqs = {}
            cov = {}
            circ = {}
            rotated = {}
            rotation_marker = {}

            self_linked_edges = set()
            single_edge_paths = {}
            missing_sequences = []

            # Parse the GFA once and collect three pieces of information:
            # S records: edge sequences and per-edge tags
            # L records: self-links, used as circularity evidence
            # P records: contig paths, used to identify single-edge contigs
            with open(input.gfa) as handle:
                for line in handle:
                    line = line.rstrip("\n")
                    if not line:
                        continue

                    fields = line.split("\t")
                    record_type = fields[0]

                    if record_type == "S":
                        edge = fields[1]
                        seq = fields[2]

                        edge_order.append(edge)

                        # Flye GFA normally embeds sequences in S records.
                        # GraphMB edge-mode needs these sequences to build assembly.fasta.
                        if seq == "*":
                            missing_sequences.append(edge)
                            continue

                        seqs[edge] = seq
                        cov[edge] = "*"
                        circ[edge] = "N"
                        rotated[edge] = "N"
                        rotation_marker[edge] = "*"

                        # dp:i:<value> is Flye's edge coverage tag.
                        # RT:z:<marker> is added by dnaapler when it rotates a sequence.
                        for tag in fields[3:]:
                            if tag.startswith("dp:i:"):
                                cov[edge] = tag.split(":", 2)[2]
                            elif tag.startswith("RT:z:"):
                                rotated[edge] = "Y"
                                rotation_marker[edge] = tag.split(":", 2)[2]

                    elif record_type == "L" and len(fields) >= 6:
                        from_edge = fields[1]
                        to_edge = fields[3]
                        overlap = fields[5]

                        # A self-link is one part of the conservative circularity call.
                        # We combine this with a single-edge P path below.
                        if from_edge == to_edge and overlap == "0M":
                            self_linked_edges.add(from_edge)

                    elif record_type == "P" and len(fields) >= 3:
                        path_name = fields[1]
                        path = fields[2]
                        path_edges = path.split(",")

                        # Flye reports simple circular contigs like:
                        # P contig_1250 edge_1250+ *
                        # These are easy to map back to one graph edge.
                        if len(path_edges) == 1:
                            edge = path_edges[0].rstrip("+-")
                            single_edge_paths[edge] = path_name

            if not edge_order:
                raise ValueError("No GFA S records found")

            if missing_sequences:
                examples = ", ".join(missing_sequences[:20])
                raise ValueError(
                    "GFA contains S records without embedded sequence. "
                    f"Examples: {examples}"
                )

            # Mark circularity after parsing all records.
            # We require both: a self-link and a single-edge contig path.
            # This avoids marking arbitrary self-linked graph components as circular contigs.
            for edge in edge_order:
                if edge in self_linked_edges and edge in single_edge_paths:
                    circ[edge] = "Y"

            with open(output.fasta, "w") as out_fa, open(output.info, "w") as out_info:
                out_info.write(
                    "#seq_name\tlength\tcov.\tcirc.\t"
                    "rotated_by_dnaapler\trotation_marker\n"
                )

                for edge in edge_order:
                    seq = seqs[edge]

                    # Write rotation information in fasta description
                    header = edge
                    if rotated[edge] == "Y":
                        header += f" rotated=True rotated_gene={rotation_marker[edge]}"
                    out_fa.write(f">{header}\n")
                    out_fa.write("\n".join(textwrap.wrap(seq, 80)) + "\n")

                    out_info.write(
                        f"{edge}\t"
                        f"{len(seq)}\t"
                        f"{cov[edge]}\t"
                        f"{circ[edge]}\t"
                        f"{rotated[edge]}\t"
                        f"{rotation_marker[edge]}\n"
                    )

            with open(log[0], "w") as lg:
                lg.write(f"Input GFA: {input.gfa}\n")
                lg.write(f"Output FASTA: {output.fasta}\n")
                lg.write(f"Output graph info: {output.info}\n")
                lg.write(f"S records: {len(edge_order)}\n")
                lg.write(f"Self-linked edges: {len(self_linked_edges)}\n")
                lg.write(f"Single-edge paths: {len(single_edge_paths)}\n")
                lg.write(f"Circular single-edge records: {sum(v == 'Y' for v in circ.values())}\n")
                lg.write(f"dnaapler-rotated records: {sum(v == 'Y' for v in rotated.values())}\n")
                lg.write("Done\n")

    ###############################################################################################
    # Polish dnaapler-reoriented flye contigs using medaka
    ###############################################################################################

    # mem_usage
    # Regression: max_mem_kb = 1.214e06 + 4.863e-02*len(input.contigs)
    # On top, it gets a baseline 5GB and every new attempt gets 5GB more.
    # We first get the contig fasta locally so that medaka's internal logic
    # to create minimap index will stay in the shadow directory only.
    rule medaka_consensus_nodes:
        input:
            contigs = rules.dnaapler_reorient_nodes.output.fasta,
            reads   = "{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz",
        output:
            fasta   = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly.dnaapler.medaka.fasta"
        shadow:
            "minimal"
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/medaka_consensus_nodes.log"
        params:
            model = MEDAKA_INFERENCE_MODEL if MEDAKA_INFERENCE_MODEL is not None else ""
        threads:
            4
        resources:
            mem=lambda wildcards, input, attempt: 1.22 + 5e-8*len(input.contigs) + 5*attempt,
            gpu=1
        conda:
            minto_dir + "/envs/ont_polishing.yml"
        shell:
            """
            time (
                rsync {input.contigs} assembly.fasta
                medaka_consensus -i {input.reads} -d assembly.fasta -o out -t {threads} {params.model}
                rsync -a out/consensus.fasta {output.fasta}
            ) >& {log}
            """

    ###############################################################################################
    # Polish dnaapler-reoriented flye gfa edge_fasta using medaka
    ###############################################################################################

    use rule medaka_consensus_nodes as medaka_consensus_edges with:
        input:
            contigs = rules.extract_gfa_edge_info_fasta.output.fasta,
            reads   = "{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz",
        output:
            fasta   = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_graph.dnaapler.medaka.fasta"
        log:
            "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/medaka_consensus_edges.log"

    ###############################################################################################
    # Common function to mark circular contigs in fasta header based on graph info txt files
    ###############################################################################################

    def mark_circular_seq(wildcards, input, output, seq_type):
        # Step 1: Read circularity info file to find structurally circular sequences
        circular = {}
        with open(input.circ_info, "r") as f:
            for line in f:
                words = line.strip().split('\t')
                if len(words) > 3 and words[3] == "Y":
                    circular[words[0]] = 1

        # Step 2: Read Dnaapler's intermediate FASTA to harvest metadata notes
        dnaapler_notes = {}
        with open(input.dnaapler_fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Extract the base ID and any trailing description fields
                    header_parts = line.strip().replace(">", "").split(maxsplit=1)
                    base_id = header_parts[0]

                    # Save the description if Dnaapler appended one (e.g., 'rotated=True rotated_gene=terL')
                    if len(header_parts) > 1:
                        dnaapler_notes[base_id] = header_parts[1]

        # Step 3: Parse polished Medaka file, combining and reapplying metadata tags
        with open(output.final_fasta, 'w') as out_fa, open(output.name_mapping, 'w') as out_mapping:
            fiter = fasta_iter(input.medaka_fasta)
            for entry in fiter:
                header, seq = entry

                # Medaka outputs the bare ID. Ensure we match the original baseline string
                base_id = header.split()[0]

                # Convert 'edge/contig' to 'EDGE/NODE' to maintain MIntO's standard layout naming scheme
                if seq_type == "EDGE":
                    seq_id = base_id.replace("edge", "EDGE")
                else:
                    seq_id = base_id.replace("contig", "NODE")
                new_header = f"{seq_id}_length_{len(seq)}"

                # A. Re-apply the structural Flye circularity flag if true
                if base_id in circular:
                    new_header = f"{new_header}_circularA"

                # B. Re-append dnaapler's functional notes if they exist
                if base_id in dnaapler_notes:
                    new_header = f"{new_header} {dnaapler_notes[base_id]}"

                # Final header!
                new_header = f"MetaFlye.{wildcards.assembly_preset}.{wildcards.nanopore}_{new_header}"

                # Write mapping
                out_mapping.write(f"{base_id}\t{new_header}\n")

                # Write sequence out
                out_fa.write(f">{new_header}\n")
                out_fa.write(seq+"\n")

    ###############################################################################################
    # This rule identifies contigs marked circular in Flye-generated assembly_info.txt.
    # It appends '_circularA' to the fasta header.
    # For dnaapler rotated contigs, it also retains the fasta comment added by dnaapler
    ###############################################################################################
    rule mark_circular_flye_nodes:
        localrule: True
        input:
            medaka_fasta   = rules.medaka_consensus_nodes.output.fasta,
            dnaapler_fasta = rules.dnaapler_reorient_nodes.output.fasta,
            circ_info      = rules.nanopore_assembly_metaflye.output.info,
        output:
            final_fasta    = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}.assembly.nodes.fasta",
            name_mapping   = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}.assembly.nodes.renaming.tsv",
        run:
            mark_circular_seq(wildcards, input, output, seq_type="NODE")

    ###############################################################################################
    # This rule identifies edges marked circular in MIntO-generated assembly_graph_info.txt.
    # It appends '_circularA' to the fasta header.
    # For dnaapler rotated contigs, it also retains the fasta comment added by dnaapler
    ###############################################################################################
    rule mark_circular_flye_edges:
        localrule: True
        input:
            medaka_fasta   = rules.medaka_consensus_edges.output.fasta,
            dnaapler_fasta = rules.extract_gfa_edge_info_fasta.output.fasta,
            circ_info      = rules.extract_gfa_edge_info_fasta.output.info,
        output:
            final_fasta    = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}.assembly.edges.fasta",
            name_mapping   = "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}.assembly.edges.renaming.tsv",
        run:
            mark_circular_seq(wildcards, input, output, seq_type="EDGE")
