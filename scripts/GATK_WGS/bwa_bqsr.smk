
rule bqsr_apply:
    input:
        ref = config["Reference"],
        bam = lambda wildcards: ("{out}/{bam}/{sample}_sorted.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample=wildcards.sample)),
        grp = lambda wildcards: ("{out}/{tmp}/{sample}_bqsr.grp".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample=wildcards.sample)),
        gatk = config["GATK"]
    output:
        bam = "{out}/{bam}/{sample}_recal.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}"),
        bai = "{out}/{bam}/{sample}_recal.bai".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}")
    threads:
        4
    shell:
        """
        java -jar {input.gatk} -T PrintReads \
            -nct {threads} \
            -R {input.ref} \
            -I {input.bam} \
            -BQSR {input.grp} \
            -o {output.bam};
        sleep 5
        """

rule bqsr_learn:
    input:
        ref = config["Reference"],
        bam = lambda wildcards: ("{out}/{bam}/{sample}_sorted.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample=wildcards.sample)),
        bai = lambda wildcards: ("{out}/{bam}/{sample}_sorted.bam.bai".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample=wildcards.sample)),
        gatk = config["GATK"]
    output:
        temp("{out}/{tmp}/{sample}_bqsr.grp".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}"))
    threads:
        4
    params:
        dbsnp = config["resource"]["dbsnp"],
        mills = config["resource"]["mills"],
        onekg = config["resource"]["onekgenomes"]
    shell:
        """
        java -jar {input.gatk} -T BaseRecalibrator \
            -nct {threads} \
            -R {input.ref} \
            -I {input.bam} \
            -knownSites {params.dbsnp} \
            -knownSites {params.mills} \
            -knownSites {params.onekg} \
            -o {output}
        """

rule index:
    input:
        lambda wildcards: "{out}/{bam}/{sample}_sorted.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample=wildcards.sample)
    output:
        "{out}/{bam}/{sample}_sorted.bam.bai".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}")
    threads:
        4
    shell:
        "sambamba index -t {threads} {input[0]}; "
        "sleep 2"

rule merge_bam:
    input:
        all = lambda wildcards : ("{out}/{tmp}/{sample}_{idx}_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample=wildcards.sample, idx=idx) for idx in range(len(config["Samples"][wildcards.sample]))),
        split = lambda wildcards : ("{out}/{tmp}/{sample}_{idx}_split_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample=wildcards.sample, idx=idx) for idx in range(len(config["Samples"][wildcards.sample]))),
        disc = lambda wildcards : ("{out}/{tmp}/{sample}_{idx}_disc_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample=wildcards.sample, idx=idx) for idx in range(len(config["Samples"][wildcards.sample])))
    output:
        all = "{out}/{bam}/{sample}_sorted.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}"),
        split = "{out}/{bam}/{sample}_split_sorted.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}"),
        disc = "{out}/{bam}/{sample}_disc_sorted.bam".format(out=config["dir"]["out"], bam=config["dir"]["bam"], sample="{sample}")
    shell:
        "samtools merge {output.all} {input.all}; "
        "samtools merge {output.disc} {input.disc}; "
        "samtools merge {output.split} {input.split}; "
        "sleep 5"

rule bwa:
    input:
        config["Reference"],
        lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["fq1"],
        lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["fq2"]
    params:
        sample = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["sample"],
        machine = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["machine"],
        flowcell = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["flowcell"],
        lane = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["lane"],
        library = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["library"],
        technology = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["technology"],
        tmp_dir=config["dir"]["out"] + '/' + config["dir"]["tmp"],
        sbx = lambda wildcards : config["Samples"][wildcards.sample][int(wildcards.idx)]["samp_barcode"]
    output:
        "{out}/{tmp}/{sample}_{idx}_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}", idx="{idx}"),
        "{out}/{tmp}/{sample}_{idx}_split_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}", idx="{idx}"),
        "{out}/{tmp}/{sample}_{idx}_disc_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}", idx="{idx}")
    threads:
        16
    run:
        read_group = r"@RG\tID:{}\tSM:{}\tLB:{}\tPL:{}\tPU:{}"
        read_group = read_group.format(
                params.flowcell + '.' + params.lane,
                params.sample,
                params.library,
                params.technology,
                params.flowcell + '.' + params.lane + '.' + params.sbx)
        print(r"""bwa mem -Y -t {threads} -R '{read_group}' {input[0]} {input[1]} {input[2]} | \
        samblaster --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
        --splitterFile \
        >(gawk '{{ if ($0~"^@") {{ print; next}} else {{ $10="*"; $11="*"; print }} }}' OFS="\t" - | \
        sambamba view -S -f bam -l 0 /dev/stdin | \
        sambamba sort -t 4 -m 1G --tmpdir={params.tmp_dir} -o {output[1]} /dev/stdin) \
        --discordantFile \
        >(gawk '{{ if ($0~"^@") {{ print; next}} else {{ $10="*"; $11="*"; print }} }}' OFS="\t" - | \
        sambamba view -S -f bam -l 0 /dev/stdin | \
        sambamba sort -t 4 -m 1G --tmpdir={params.tmp_dir} -o {output[2]} /dev/stdin) | \
        sambamba view -S -f bam -l 0 /dev/stdin | \
        sambamba sort -t {threads} -m 60G --tmpdir={params.tmp_dir} -o {output[0]} /dev/stdin; \
        sleep 10""")

        shell(
        r"""bwa mem -Y -t {threads} -R '{read_group}' {input[0]} {input[1]} {input[2]} | \
        samblaster --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
        --splitterFile \
        >(gawk '{{ if ($0~"^@") {{ print; next}} else {{ $10="*"; $11="*"; print }} }}' OFS="\t" - | \
        sambamba view -S -f bam -l 0 /dev/stdin | \
        sambamba sort -t 4 -m 1G --tmpdir={params.tmp_dir} -o {output[1]} /dev/stdin) \
        --discordantFile \
        >(gawk '{{ if ($0~"^@") {{ print; next}} else {{ $10="*"; $11="*"; print }} }}' OFS="\t" - | \
        sambamba view -S -f bam -l 0 /dev/stdin | \
        sambamba sort -t 4 -m 1G --tmpdir={params.tmp_dir} -o {output[2]} /dev/stdin) | \
        sambamba view -S -f bam -l 0 /dev/stdin | \
        sambamba sort -t {threads} -m 60G --tmpdir={params.tmp_dir} -o {output[0]} /dev/stdin; \
        sleep 10""")

