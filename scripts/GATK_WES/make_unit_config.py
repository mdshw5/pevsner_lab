#!/usr/bin/env python3

import sys
version = sys.version_info
if version[0] < 3 or version[1] < 3:
    sys.exit("This script requires at least python 3.3")

import argparse
import json
import string
import random
import os.path
import shutil
import subprocess
import os
import re

def test_env(commands):
    for cmd in commands:
        if not shutil.which(cmd):
            sys.exit("Error: {cmd} not found in the PATH".format(cmd=cmd))

def try_python2_import(modules):
    for mod in modules:
        try:
            subprocess.check_call("python2 -c \"import imp; imp.find_module('{mod}')\"".format(mod=mod), shell=True)
        except subprocess.CalledProcessError:
            sys.exit("Error: requires python2 with the {mod} module".format(mod=mod))

def genomestrip_test_env():
    test_env(["Rscript", "bwa", "samtools", "tabix"])
    sv_dir = os.environ.get("SV_DIR")
    if not sv_dir:
        sys.exit("Error: please set SV_DIR to the base directory of the GenomeStrip installation")
    classpath = os.environ.get("GENOSTRIP_CLASSPATH")
    if not classpath:
        sys.exit("Error: please set GENOSTRIP_CLASSPATH before configuration")

def manta_test_env():
    test_env(["samtools", "bwa", "samblaster", "sambamba", "configManta.py"])
    try:
        subprocess.check_call("""python -c $'import sys\nif sys.version_info[0] > 2:   sys.exit(1)'""", shell=True)
    except subprocess.CalledProcessError:
        sys.exit("Error: the default python executable must be python2")

def lumpy_test_env():
    test_env(["samtools", "bwa", "samblaster", "sambamba", "lumpy", "lumpyexpress"])
    try_python2_import(["scipy", "numpy", "pysam"])

def gatk_test_env():
    test_env(["samtools", "bwa", "samblaster", "sambamba", "bcftools"])

def fermikit_test_env():
    test_env(["samtools", "bfc", "seqtk", "gzip", "ropebwt2", "fermi2", "bwa", "sambamba", "htsbox", "k8"])

def freebayes_test_env():
    test_env(["samtools", "bwa", "samblaster", "sambamba", "freebayes", "bgzip", "tabix"])

def parse_args():
    parser = argparse.ArgumentParser(description="Create a json config file for snakemake",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--fastq", type=str, nargs='+', help="The input fastq files (read 1)")
    parser.add_argument("--reference", type=str, default="/mnt/data/reference/hs37d5.fa", help="The reference aligned to.")
    parser.add_argument("--technology", type=str, default="ILLUMINA", help="The technology used to produce the data")
    parser.add_argument("--out_dir", type=str, default="./Data_Files", help="The output directory for the files")
    parser.add_argument("--config", help="The config file to write to")
    subparsers = parser.add_subparsers(dest="subparser_name")

    # The GATK parser #
    gatk_parser = subparsers.add_parser("GATK")
    gatk_parser.add_argument("--n_hc_intervals", type=int, default=48, help="The number of intervals to use when running the haplotype caller in parallel.")
    gatk_parser.add_argument("--n_geno_intervals", type=int, default=48, help="The number of intervals to use when running GenotypeGVCFs")
    gatk_parser.add_argument("--ploidy", type=int, default=2, help="The sample ploidy")
    gatk_parser.add_argument("--bam_dir", type=str, default="BAM", help="The output subdirectory for the processed BAM files")
    gatk_parser.add_argument("--gvcf_dir", type=str, default="GVCF", help="The output subdirectory for the GVCF files")
    gatk_parser.add_argument("--vcf_dir", type=str, default="VCF", help="The output subdirectory for the VCF files")
    gatk_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for sambamba sorting and file subsets temporary files")
    gatk_parser.add_argument("--dbsnp", type=str, default="/mnt/data/reference/GATK_b37/dbsnp_138.b37.vcf", help="A file of dbSNP variants to use as a reference.")
    gatk_parser.add_argument("--hapmap", type=str, default="/mnt/data/reference/GATK_b37/hapmap_3.3.b37.vcf", help="A file of hapmap SNPs to use in recalibration.")
    gatk_parser.add_argument("--omni", type=str, default="/mnt/data/reference/GATK_b37/1000G_omni2.5.b37.vcf", help="A file of OMNI SNPs to use in recalibration.")
    gatk_parser.add_argument("--onekgenomes", type=str, default="/mnt/data/reference/GATK_b37/1000G_phase1.snps.high_confidence.b37.vcf", help="A file of 1000 genomes SNPs to use in recalibration.")
    gatk_parser.add_argument("--mills", type=str, default="/mnt/data/reference/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf", help="A file of Mills et al. validated indels to use in recalibration.")
    gatk_parser.add_argument("gatk", type=str, help="The GATK .jar")
    gatk_parser.add_argument("rm_dup_script", type=str, help="A script to remove duplicated variants at interval edges")

    # The fermikit parser #
    fermikit_parser = subparsers.add_parser("fermikit")
    fermikit_parser.add_argument("--bam_dir", type=str, default="BAM", help="The output subdirectory for the processed BAM files")
    fermikit_parser.add_argument("--vcf_dir", type=str, default="VCF", help="The output subdirectory for the VCF files")
    fermikit_parser.add_argument("--ec_dir", type=str, default="EC", help="The output subdirectory for error-corrected fastq files.")
    fermikit_parser.add_argument("--asmbl_dir", type=str, default="Asmbl", help="The output subdirectory for the de novo assembly")
    fermikit_parser.add_argument("--sv_dir", type=str, default="SV", help="The output subdirectory for detected structural variants")
    fermikit_parser.add_argument("--tmp_dir", type=str, default="tmp", help="The output subdirectory for temporary files")
    fermikit_parser.add_argument("hapdip", help="The hapdip javascript program")
    fermikit_parser.add_argument("read_length", type=int, help="The length of the sequence reads")

    # The freebayes parser #
    freebayes_parser = subparsers.add_parser("freebayes")
    freebayes_parser.add_argument("--n_freebayes_intervals", type=int, default=12, help="The number of intervals to use when running freebayes in parallel")
    freebayes_parser.add_argument("--bam_dir", type=str, default="BAM", help="The output subdirectory for the processed BAM files")
    freebayes_parser.add_argument("--vcf_dir", type=str, default="VCF", help="The output subdirectory for the VCF files")
    freebayes_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for sambamba sorting and file subsets temporary files")
    freebayes_parser.add_argument("rm_dup_script", type=str, help="A script to remove duplicated variants at interval edges")
    # For VQSR
    freebayes_parser.add_argument('--dbsnp', type=str, default='/mnt/data/reference/GATK_b37/dbsnp_138.b37.vcf', help='dbSNP variants for recalibration')
    freebayes_parser.add_argument('--hapmap', type=str, default='/mnt/data/reference/GATK_b37/hapmap_3.3.b37.vcf', help='Hapmap SNPs for recalibration')
    freebayes_parser.add_argument("--omni", type=str, default="/mnt/data/reference/GATK_b37/1000G_omni2.5.b37.vcf", help="A file of OMNI SNPs to use in recalibration.")
    freebayes_parser.add_argument("--onekgenomes", type=str, default="/mnt/data/reference/GATK_b37/1000G_phase1.snps.high_confidence.b37.vcf", help="A file of 1000 genomes SNPs to use in recalibration.")
    freebayes_parser.add_argument("--mills", type=str, default="/mnt/data/reference/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf", help="A file of Mills et al. validated indels to use in recalibration.")
    freebayes_parser.add_argument("--n_hc_intervals", type=int, default=48, help="The number of intervals to use when running the haplotype caller in parallel.")
    freebayes_parser.add_argument("--n_geno_intervals", type=int, default=48, help="The number of intervals to use when running GenotypeGVCFs")
    freebayes_parser.add_argument('gatk', type=str, help='GATK .jar')

    # The lumpy parser #
    lumpy_parser = subparsers.add_parser("lumpy")
    lumpy_parser.add_argument("--bam_dir", type=str, default="BAM", help="The output subdirectory for the processed BAM files")
    lumpy_parser.add_argument("--vcf_dir", type=str, default="VCF", help="The output subdirectory for the VCF files")
    lumpy_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for temporary files")
    lumpy_parser.add_argument("svtyper", type=str, help="The svtyper script")
    lumpy_parser.add_argument("lumpy_exclude", type=str, help="A script to create lumpy exclude regions")

    # The manta parser #
    manta_parser = subparsers.add_parser("manta")
    manta_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for temporary files")
    manta_parser.add_argument("--bam_dir", type=str, default="BAM", help="The output subdirectory for the processed BAM files")
    manta_parser.add_argument("--manta_dir", type=str, default="Manta", help="The output subdirectory for manta")

    # The GenomeSTRIP parser #
    genomestrip_parser = subparsers.add_parser("GenomeSTRIP")
    genomestrip_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for temporary files")
    genomestrip_parser.add_argument("--bam_dir", type=str, default="BAM", help="The output subdirectory for the processed BAM files")
    genomestrip_parser.add_argument("--cnv_dir", type=str, default="CNV", help="A subdirectory for GenomeSTRIP output")
    genomestrip_parser.add_argument("--metadata_dir", type=str, default="Strip_Meta", help="A subdirectory for GenomeSTRIP sample metadata")
    genomestrip_parser.add_argument("--log_dir", type=str, default="Strip_Logs", help="A subdirectory for GenomeSTRIP log files")
    genomestrip_parser.add_argument("--coverage", type=int, default=30, help="Average coverage of all samples (only useful for low coverage genomes")
    genomestrip_parser.add_argument("check_sex", type=str, help="A script for determining sample sex")
    genomestrip_parser.add_argument("interval_creator", type=str, help="A script for creating a list of chromosomes")

    return parser.parse_args()

def main_genomestrip(args):
    genomestrip_test_env()

    args.directories.update( {
            "tmp": args.tmp_dir,
            "bam": args.bam_dir,
            "cnv": args.cnv_dir,
            "log": args.log_dir,
            "metadata": args.metadata_dir } )

    tiling_size = 1000
    tiling_overlap = 500
    reference_gap_length = 1000
    boundary_precision = 100
    min_refined_length = 500

    if args.coverage < 30:
        tiling_size = int(tiling_size * 30 / args.coverage)
        tiling_overlap = int(tiling_overlap * 30 / args.coverage)
        reference_gap_length = int(reference_gap_length * 30 / args.coverage / 2)
        boundary_precision = int(boundary_precision * 30 / args.coverage / 2.5)
        min_refined_length = int(min_refined_length * 30 / args.coverage)

    sv_dir = os.environ.get("SV_DIR")
    classpath = os.environ.get("GENOSTRIP_CLASSPATH")

    args.configs.update( {
            "size": tiling_size,
            "overlap": tiling_overlap,
            "gap": reference_gap_length,
            "boundary": boundary_precision,
            "refined": min_refined_length,
            "check_sex": os.path.abspath(args.check_sex),
            "interval_creator": os.path.abspath(args.interval_creator),
            "classpath": classpath,
            "sv_dir": sv_dir,
            "dir": args.directories } )

def main_manta(args):
    manta_test_env()

    args.directories.update( {
            "tmp": args.tmp_dir,
            "bam": args.bam_dir,
            "manta": args.manta_dir } )

    args.configs.update( {
            "dir": args.directories } )

def main_lumpy(args):
    lumpy_test_env()

    args.directories.update( {
            "tmp": args.tmp_dir,
            "bam": args.bam_dir,
            "vcf": args.vcf_dir } )

    args.configs.update( {
            "dir": args.directories,
            "svtyper": args.svtyper,
            "lumpy_exclude": args.lumpy_exclude } )
        
def main_gatk(args):
    gatk_test_env()
    
    args.directories.update( {
            "tmp": args.tmp_dir,
            "bam": args.bam_dir,
            "gvcf": args.gvcf_dir,
            "vcf": args.vcf_dir } )

    resource = {
        "dbsnp": args.dbsnp,
        "hapmap":args.hapmap,
        "omni":args.omni,
        "onekgenomes":args.onekgenomes,
        "mills":args.mills}

    args.configs.update( {
            "GATK": args.gatk,
            "n_hc_intervals": str(args.n_hc_intervals),
            "n_geno_intervals": str(args.n_geno_intervals),
            "resource": resource,
            "ploidy": str(args.ploidy),
            "dir": args.directories,
            "rmdup": os.path.abspath(args.rm_dup_script)} )

def main_fermikit(args):
    fermikit_test_env()
    
    min_unitig_ovlp = str(int(51 + (args.read_length - 101) * 0.4 + 0.499))
    min_merge_ovlp = str(int(args.read_length * 0.75) + 1)
    min_clean_ovlp = str(int(min_unitig_ovlp) + 5)

    args.directories.update( {
            "ec": args.ec_dir,
            "asmbl": args.asmbl_dir,
            "bam": args.bam_dir,
            "vcf": args.vcf_dir,
            "sv": args.sv_dir,
            "tmp": args.tmp_dir} )

    args.configs.update( {
            "dir": args.directories,
            "hapdip": args.hapdip,
            "min_unitig_ovlp": min_unitig_ovlp,
            "min_merge_ovlp": min_merge_ovlp,
            "min_clean_ovlp": min_clean_ovlp} )

def main_freebayes(args):
    freebayes_test_env()

    args.directories.update( {
            "bam": args.bam_dir,
            "vcf": args.vcf_dir,
            "tmp": args.tmp_dir } )

    resource = {
            'dbsnp': args.dbsnp,
            'hapmap': args.hapmap,
            'omni': args.omni,
            'onekgenomes': args.onekgenomes,
            'mills': args.mills}

    args.configs.update( {
            "dir": args.directories,
            'resource': resource,
            'GATK': args.gatk,
            "n_freebayes_intervals": args.n_freebayes_intervals,
            'n_hc_intervals': str(args.n_hc_intervals),
            'n_geno_intervals': str(args.n_geno_intervals),
            "rmdup": os.path.abspath(args.rm_dup_script)} )

if __name__ == "__main__":
    args = parse_args()

    samples = {}
    for i, fq1 in enumerate(args.fastq):
        if fq1.endswith("R1.fastq.gz"):
            fq2 = fq1[:-11] + "R2.fastq.gz"
        else:
            sys.exit("Can not find fastq2")

        fq1 = os.path.abspath(fq1)
        fq2 = os.path.abspath(fq2)
        if not os.path.isfile(fq1):
            sys.exit("{} is not a file".format(fq1))
        if not os.path.isfile(fq2):
            sys.exit("{} is not a file".format(fq2))

        m = re.match(r".*/([\w-]+)_([\da-zA-Z-]+)_([\dA-Za-z-]+)_([\dA-Za-z]+)_([ATCGN]+)_([\d]+)_R1.fastq.gz", fq1) # sample, library, machine, flowcell, sample_barcode, lane
        sample = m.group(1)
        library = m.group(2)
        machine = m.group(3)
        flowcell = m.group(4)
        sbx = m.group(5)
        lane = m.group(6)
        
        if not sample in samples:
            samples[sample] = []

        cur_sample = {
            "sample": sample,
            "fq1": fq1,
            "fq2": fq2,
            "library": library,
            "machine": machine,
            "flowcell": flowcell,
			"samp_barcode": sbx,
            "lane": lane,
            "technology": args.technology }
        samples[sample].append(cur_sample)

    args.directories = {"out": os.path.abspath(args.out_dir)}

    args.configs = { "Reference": args.reference,
                     "Samples": samples}

    if args.subparser_name == "GATK":
        main_gatk(args)
    elif args.subparser_name == "fermikit":
        main_fermikit(args)
    elif args.subparser_name == "freebayes":
        main_freebayes(args)
    elif args.subparser_name == "lumpy":
        main_lumpy(args)
    elif args.subparser_name == "manta":
        main_manta(args)
    elif args.subparser_name == "GenomeSTRIP":
        main_genomestrip(args)
    else:
        sys.exit("Unrecognized command")
    
    with open(args.config, 'w') as f:
        json.dump(args.configs, f, indent=2)
