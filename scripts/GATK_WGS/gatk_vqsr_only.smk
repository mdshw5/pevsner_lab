#v1
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
        vcf="{out}/{vcf}/combined_filt.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"])

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
        "{out}/{vcf}/combined_indels.tranches".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
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
        "{out}/{vcf}/combined_indels.tranches".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
    threads:
        6
    shell:
        "java -Xmx20g -jar {input.gatk} -T VariantRecalibrator -nt {threads} -R {input.ref} -input {input[0]} "
        "-resource:mills,known=true,training=true,truth=true,prior=12.0 {input.mills} "
        "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} "
        "--mode INDEL --recal_file {output[0]} --tranches_file {output[1]} -an QD -an FS -an ReadPosRankSum -an SOR -an MQRankSum "
        "-an MQ -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 -tranche 99.95 -tranche 99.93 -tranche 99.90 "
        "-tranche 99.8 -tranche 99.5 -tranche 99 -tranche 98 -tranche 97.5 -tranche 97 -tranche 96 -tranche 95 -tranche 90 -mG 5 "

rule snp_apply:
    input:
        "{out}/{vcf}/combined_raw.vcf".format(out=config["dir"]["out"], vcf=config["dir"]["vcf"]),
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
        "-an QD -an FS -an ReadPosRankSum -an SOR -an MQRankSum -an MQ -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 "
        "-tranche 99.95 -tranche 99.93 -tranche 99.90 -tranche 99.8 -tranche 99.5 -tranche 99 -tranche 98 -tranche 90 -mG 5"
