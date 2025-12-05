#!/usr/bin/env python

import os

include: 'resources.smk'

################################################################################################
# Wrapper to create mapper-specific index files for any given fasta file, with extensions .fna or .fasta.
################################################################################################

def get_fasta_index_path(fasta, mapper):
    # Get UC mapper name
    my_mapper = mapper.upper()

    # Split the fasta file name
    basedir, filename = os.path.split(fasta)

    # Enumerate the index files
    if my_mapper == "BWA":
        ext_list = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
    elif my_mapper == "STROBEALIGN":
        ext_list = ['']
    else:
        raise Exception(f"MIntO error: Unexpected value: MAPPER='{my_mapper}'. Must be one of {BWA, STROBEALIGN}")

    # Return a list of the index files
    index_files = expand("{location}/{mapname}_index/{filename}{ext}",
                         location = basedir,
                         mapname  = my_mapper,
                         filename = filename,
                         ext      = ext_list)
    return(index_files)

############################################
# bwa-mem2 index can handle gzipped files
############################################

# Memory requirements:
# --------------------
# Baseline    : 5 GB
# Size-based  : 22 byte per byte file size (3x22=66 if gzipped)
# New attempts: +40 GB each time
rule BWA_index:
    input:
        fasta="{somewhere}/{something}.{fasta}"
    output:
        "{somewhere}/BWA_index/{something}.{fasta}.0123",
        "{somewhere}/BWA_index/{something}.{fasta}.amb",
        "{somewhere}/BWA_index/{something}.{fasta}.ann",
        "{somewhere}/BWA_index/{something}.{fasta}.bwt.2bit.64",
        "{somewhere}/BWA_index/{something}.{fasta}.pac"
    log:
        "{somewhere}/BWA_index/BWA_index.{something}.{fasta}.log"
    wildcard_constraints:
        fasta     = r'fasta|fna|fasta\.gz|fna\.gz',
        something = r'[^/]+'
    shadow:
        "minimal"
    resources:
        mem = lambda wildcards, input, attempt: 5 + int((66 if input.fasta.endswith('.gz') else 22)*get_file_size_gb(input.fasta)) + 40*(attempt-1),
    threads: 4
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """
        outdir=$(dirname {output[0]})
        prefix=$(basename {input.fasta})
        time (
            bwa-mem2 index {input.fasta} -p $prefix
            rsync -av --itemize-changes $prefix.* $outdir/
        ) >& {log}
        """

############################################
# strobealign index can handle gzipped files
############################################

# Memory requirements:
# --------------------
# Baseline    : 5 GB
# Size-based  : 13 byte per byte file size gzipped + 20
# New attempts: +40 GB each time
rule STROBEALIGN_index:
    input:
        fasta="{somewhere}/{something}.{fasta}",
        meanlen_txt=f"{working_dir}/output/2-qc/{omics}.mean_length.txt"
    output:
        fasta="{somewhere}/STROBEALIGN_index/{something}.{fasta}"
    log:
        "{somewhere}/STROBEALIGN_index/sba_index.{something}.{fasta}.log"
    wildcard_constraints:
        fasta     = r'fasta|fna|fasta\.gz|fna\.gz',
        something = r'[^/]+'
    resources:
        mem = lambda wildcards, input, attempt: 5 + int((13 if input.fasta.endswith('.gz') else 4)*get_file_size_gb(input.fasta)) + 40*(attempt-1),
    threads: 8
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #strobealign
    shell:
        """
        time (
            r_arg="$(cat {input.meanlen_txt})"
            # strobealign does realpath on input file and places sti there.
            # so hardlink source fasta file into the same directory as index, 
            # since softlink will resolve to original location.
            ln --force {input.fasta} {output.fasta}
            strobealign --create-index -t {threads} {output.fasta} -r $r_arg
        ) >& {log}
        """
