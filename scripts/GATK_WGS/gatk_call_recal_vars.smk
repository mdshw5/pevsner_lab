include: "bwa_bqsr.smk"
include: "make_directories.smk"

def get_intervals(ref_file, n_intervals):
    import re
    import math

    n_intervals = int(n_intervals)
    # Calculate step size #
    total_size = 0
    chroms = []
    bp = []

    with open(ref_file + '.fai') as f:
        for line in f:
            line = line.rstrip().split()
            total_size += int(line[1])
            chroms.append(line[0])
            bp.append(int(line[1]))
    step_size = math.ceil(total_size / n_intervals)

    # Get intervals #
    intervals = {}
    start_pos = []
    last_chrom = ''
    last_pos = 1
    passed = 0
    cur_interval = []

    for i in range(len(chroms)):
        if chroms[i] == 'hs37d5': # Skip
            continue
        while (last_pos + step_size - passed < bp[i]):
            cur_interval.append('-L ' + chroms[i] + ':' + str(last_pos) + '-' + str(last_pos + step_size - passed))
            start_pos.append(re.split(r"[:-]", cur_interval[0][3:])[:2])
            start_pos[-1] = '_'.join(start_pos[-1])
            intervals[start_pos[-1]] = ' '.join(cur_interval)
            cur_interval = []
            last_pos += step_size - passed
            passed = 0
        else:
            cur_interval.append('-L ' + chroms[i] + ':' + str(last_pos) + '-' + str(bp[i]))
            passed += bp[i] - last_pos
            last_pos = 1
    if cur_interval:
        start_pos.append(re.split(r"[:-]", cur_interval[0][3:])[:2])
        start_pos[-1] = '_'.join(start_pos[-1])
        intervals[start_pos[-1]] = ' '.join(cur_interval)

    return (intervals, start_pos)

intervals, start_pos = get_intervals(config["Reference"], config["n_hc_intervals"])
geno_intervals, geno_start = get_intervals(config["Reference"], config["n_geno_intervals"])

rule all:
    input:
        vcf="{out}/{vcf}/combined_raw.vcf.gz".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])

rule filter:
    input:
        "{out}/{vcf}/combined_recal.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])
    output:
        "{out}/{vcf}/combined_filt.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])
    shell:
        r"grep '^#\|PASS' {input} > {output}"

rule indel_apply:
    input:
        "{out}/{vcf}/combined_snp_recal.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        "{out}/{vcf}/combined_indels.recal".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        "{out}/{vcf}/combined_indles.tranches".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        gatk=config['GATK'],
        ref=config['Reference']
    output:
        "{out}/{vcf}/combined_recal.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])
    threads:
        4
    shell:
        "java -Xmx20g -jar {input.gatk} -T ApplyRecalibration -nt {threads} -R {input.ref} -input {input[0]} --mode INDEL -ts_filter_level 98.0 "
        "-recalFile {input[1]} -tranchesFile {input[2]} -o {output}"

rule indel_recal:
    input:
        "{out}/{vcf}/combined_snp_recal.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        gatk=config['GATK'],
        ref=config['Reference'],
        dbsnp=config['resource']['dbsnp'],
        mills=config['resource']['mills']
    output:
        "{out}/{vcf}/combined_indels.recal".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        "{out}/{vcf}/combined_indelss.tranches".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
    threads:
        6
    shell:
        "java -Xmx20g -jar {input.gatk} -T VariantRecalibrator -nt {threads} -R {input.ref} -input {input[0]} "
        "-resource:mills,known=true,training=true,truth=true,prior=12.0 {input.mills} "
        "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} "
        "--mode INDEL --recal_file {output[0]} --tranches_file {output[1]} -an DP -an QD -an FS -an ReadPosRankSum -an SOR -an MQRankSum "
        "-an MQ -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 -tranche 99.95 -tranche 99.93 -tranche 99.90 "
        "-tranche 99.8 -tranche 99.5 -tranche 99 -tranche 98 -tranche 97.5 -tranche 97 -tranche 96 -tranche 95 -tranche 90 -mG 5 "

rule snp_apply:
    input:
        "{out}/{vcf}/combined_snps.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        "{out}/{vcf}/combined_snps.recal".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        "{out}/{vcf}/combined_snps.tranches".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        gatk=config['GATK'],
        ref=config['Reference']
    output:
        #"Data_Files/VCF/fam_snp_recal.vcf"
        "{out}/{vcf}/combined_snp_recal.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])
    threads:
        4
    shell:
        "java -Xmx20g -jar {input.gatk} -T ApplyRecalibration -nt {threads} -R {input.ref} "
        "-input {input[0]} --mode SNP --ts_filter_level 99.0 -recalFile {input[1]} "
        "-tranchesFile {input[2]} -o {output}"

rule snp_recal:
    input:
        "{out}/{vcf}/combined_raw.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        gatk=config['GATK'],
        ref=config['Reference'],
        dbsnp=config['resource']['dbsnp'],
        hapmap=config['resource']['hapmap'],
        omni=config['resource']['omni'],
        onekg=config['resource']['onekgenomes']
    output:
        "{out}/{vcf}/combined_snps.recal".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
        "{out}/{vcf}/combined_snps.tranches".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])
    threads:
        6
    shell:
        "java -Xmx20g -jar {input.gatk} -T VariantRecalibrator -nt {threads} -R {input.ref} "
        "-input {input[0]} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} "
        "-resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} "
        "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.onekg} "
        "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} "
        "--mode SNP --recal_file {output[0]} --tranches_file {output[1]} "
        "-an DP -an QD -an FS -an ReadPosRankSum -an SOR -an MQRankSum -an MQ -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 "
        "-tranche 99.95 -tranche 99.93 -tranche 99.90 -tranche 99.8 -tranche 99.5 -tranche 99 -tranche 98 -tranche 90 -mG 5"

def get_vcf_subsets(wildcards):
    infiles = []
    for start in geno_start:
        infile = "{out}/{tmp}/calls/{start_pos}_combined_raw.vcf"
        infile = infile.format(out=config["dir"]["out"], 
                               tmp=config["dir"]["tmp"], 
                               start_pos=start)
        infiles.append(infile)
    return infiles

rule concat_vcf:
    input:
        infiles = get_vcf_subsets,
        rmdup = config["rmdup"]
    output:
       "{out}/{vcf}/combined_raw.vcf.gz".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])
    threads:
        1
    run:
        vcf1 = input.infiles[0]
	shell('(grep "^#" {vcf1}; grep -h -v "^#" {input.infiles}) | {input.rmdup} | bgzip -c > {output}; sleep 2')

rule genotype:
    input:
        gvcfs=("{out}/{gvcf}/{sample}.g.vcf".format(out=config["dir"]["out"], gvcf=config["dir"]["gvcf"], sample=sample) for sample in config["Samples"].keys()),
        gatk=config['GATK'],
        ref=config['Reference'],
        dbsnp=config['resource']['dbsnp']
    output:
        "{out}/{tmp}/calls/{start_pos}_combined_raw.vcf".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], start_pos="{start_pos}")
    threads:
        1
    run:
        gvcfs = '-V ' + ' -V '.join(input.gvcfs)
        start = str(wildcards.start_pos)
        interval = geno_intervals[str(start)]
        shell("java -Xmx4g -jar {input.gatk} -T GenotypeGVCFs -nt {threads} "
        "-R {input.ref} {gvcfs} -D {input.dbsnp} "
        "-o {output} -stand_emit_conf 10 {interval} "
        "-G StandardAnnotation -A AlleleBalance -A TandemRepeatAnnotator "
        "-A ClippingRankSumTest -A GCContent -A MappingQualityZero -A SpanningDeletions -A StrandOddsRatio; sleep 2")

def get_gvcf_subsets(wildcards):
    infiles = []
    for start in start_pos:
        infile = "{out}/{tmp}/{sample}/{start_pos}_subset.g.vcf"
        infile = infile.format(out=config["dir"]["out"], 
                               tmp=config["dir"]["tmp"],
                               sample=wildcards.sample, 
                               start_pos=start)
        infiles.append(infile)
    return infiles

rule concat_gvcf:
    input:
        infiles=get_gvcf_subsets,
        ref=config['Reference'],
        gatk=config['GATK']
    output:
        "{out}/{gvcf}/{sample}.g.vcf".format(out=config["dir"]["out"], gvcf=config["dir"]["gvcf"], sample="{sample}")
    threads:
        1
    run:
        infiles = '-V ' + ' -V '.join(input.infiles)
        shell("java -Xmx4g -cp {input.gatk} org.broadinstitute.gatk.tools.CatVariants "
        "-R {input.ref} {infiles} -out {output} -assumeSorted; sleep 2")

rule call_vars:
    input:
        "{out}/{bam}/{sample}_recal.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}"),
        "{out}/{bam}/{sample}_recal.bai".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}"),
        gatk=config['GATK'],
        ref=config['Reference'],
    output:
        "{out}/{tmp}/{sample}/{start_pos}_subset.g.vcf".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}", start_pos="{start_pos}")
    threads:
        1
    run:
        start = str(wildcards.start_pos)
        interval = intervals[str(start)]
        outfile = "{out}/{tmp}/{sample}/{start}_subset.g.vcf".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{wildcards.sample}", start=start)
        shell(
                "java -Xmx4g -jar {input.gatk} -T HaplotypeCaller -R {input.ref} "
                "-I {input[0]} -o {output} "
                "{interval} --emitRefConfidence GVCF "
                "--variant_index_type LINEAR --variant_index_parameter 128000 "
                "-G StandardAnnotation -A AlleleBalance -A TandemRepeatAnnotator "
                "-A ClippingRankSumTest -A GCContent "
                "-A MappingQualityZero -A SpanningDeletions "
                "-A StrandOddsRatio; sleep 2"
            )

