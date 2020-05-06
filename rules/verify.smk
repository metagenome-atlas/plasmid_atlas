

VERIFY_DATA_URL ="https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/17904323/nbc_hmms.hmm.gz"
localrules: unzip,downlad_verification_tool
rule downlad_verification_tool:
    output:
        os.path.join(config['database_dir'],"CircularVerify/nbc_hmms.hmm")
    shell:
        "wget {VERIFY_DATA_URL} -O {output}.gz ; "
        " gunzip {output}.gz "


rule unzip:
    input:
        "plasmids/plasmids.deduplicated.fasta.gz"
    output:
        temp("plasmids/plasmids.deduplicated.fasta")
    shadow:
        "minimal"
    shell:
        "zcat {input} > {output}"


rule verify:
    input:
        database= rules.downlad_verification_tool.output,
        fasta= "plasmids/plasmids.deduplicated.fasta",
    log:
        "logs/plasmids/verify.log"
    params:
        script= os.path.join(os.path.dirname(workflow.snakefile),
                             "scrpts/viralveify/viralverify.py")
    benchmark:
        "logs/benchmark/plasmid/verify.txt"
    output:
        directory("plasmids/verification/")
    conda:
        "../envs/verify.yaml"
    shell:
        "python {params.script} "
        " -f {input.fasta}"
        " -o {output} "
        " --hmm input.database"
        " -p"
        " -t {threads}"
        " &> {log}"
