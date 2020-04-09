

rule deduplicate:
    input:
        expand("{sample}/plasmids/plasmids.fasta.gz",sample=SAMPLES)
    output:
        "plasmids/plasmids.deduplicated.fasta.gz"
    log:
        "logs/plasmids/deduplicate.txt"
    params:
        input= lambda wc, input: ','.join(input)
    threads:
        config['threads']
    shell:
        "dedupe.sh ow=t"
        " cluster=t pickbestrepresentative "
        " in={params.input}"
        " out={output.output}"
        " threads={threads}"
        " -Xmx{resources.mem}g"
        " 2> {log} "
