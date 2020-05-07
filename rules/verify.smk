import os

VERIFY_DATA_URL ="https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/17904323/nbc_hmms.hmm.gz"
localrules: unzip,downlad_verification_tool,link_output
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
    output:
        "plasmids/verification/plasmids.deduplicated_result_table.csv", # no header
        directory("plasmids/verification/Prediction_results_fasta"),
    log:
        "logs/plasmids/verify.log"
    params:
        script= os.path.join(os.path.dirname(workflow.snakefile),
                             "scripts/viralverify/viralverify.py"),
        outdir = "plasmids/verification"
    benchmark:
        "logs/benchmark/plasmid/verify.txt"

    conda:
        "../envs/verify.yaml"
    shell:
        "python {params.script} "
        " -f {input.fasta}"
        " -o {output} "
        " --hmm {input.database}"
        " -p"
        " -t {threads}"
        " &> {log}"

rules link_output:
    input:
        "plasmids/verification/Prediction_results_fasta"
    output:
        "plasmids/verification/circular_viruses.fasta",
        "plasmids/verification/circular_plasmids.fasta"
    run:

        from os.path import join

        output_dir= os.path.dirname(output[0])
        input_dir_rel = os.path.relpath(input[0],output_dir)


        for i,fasta in enumerate(["plasmids.deduplicated_virus.fasta","plasmids.deduplicated_plasmid.fasta"]):

            if os.path.exists(os.path.join(input[0],fasta)):

                os.symlink(join(input_dir_rel,fasta),join(output_dir,fasta))
            else:
                shell("touch {output_dir}/{fasta}")
