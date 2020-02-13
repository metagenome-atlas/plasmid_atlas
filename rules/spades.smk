
ASSEMBLY_FRACTIONS=['R1','R2','me']
assembly_preprocessing_steps='QC.errorcorr.merged'
assembly_params=dict()
assembly_params['spades'] = {'meta':'--meta','normal':''}

def spades_parameters(wc,input):
    if not os.path.exists("{sample}/assembly/params.txt".format(sample=wc.sample)):

        params={}

        reads = dict(zip(ASSEMBLY_FRACTIONS,input))

        if not PAIRED_END:
            params['inputs']= " -s {se} ".format(**reads)
        else:
            params['inputs']= " --pe1-1 {R1} --pe1-2 {R2} ".format(**reads)

            if 'se' in ASSEMBLY_FRACTIONS:
                params['inputs']+= "--pe1-s {se} ".format(**reads)
            if 'me' in ASSEMBLY_FRACTIONS:
                params['inputs']+= "--pe1-m {me} ".format(**reads)

        # Long reads:

        if (config['longread_type'] is not None) & (str(config['longread_type']).lower()!='none'):

            long_read_file = get_files_from_sampleTable(wc.sample,'longreads')[0]
            params['longreads'] = " --{t} {f} ".format(t=config['longread_type'],f=long_read_file)
        else:
            params['longreads'] = ""


        params['preset'] = assembly_params['spades'][config['spades_preset']]
        params['skip_error_correction'] = "--only-assembler" if config['spades_skip_BayesHammer'] else ""
        params['extra'] = config['spades_extra']


    else:

        params = {"inputs": "--restart-from last",
                  "preset":"",
                  "skip_error_correction":"",
                  "extra":"",
                  "longreads":""}

    params['outdir']= "{sample}/plasmids/spades".format(sample=wc.sample)

    return params


rule plasmid_spades:
    input:
        expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=ASSEMBLY_FRACTIONS,
            assembly_preprocessing_steps=assembly_preprocessing_steps)
    output:
        "{sample}/plasmids/spades/contigs.fasta",
        "{sample}/plasmids/spades/scaffolds.fasta"
    benchmark:
        "logs/benchmarks/plasmids/spades/{sample}.txt"
    params:
        p= lambda wc,input: spades_parameters(wc,input),
        k = config.get("spades_k", 'auto'),
    log:
        "{sample}/logs/plasmids/spades.log"
    conda:
        "../envs/assembly.yaml"
    threads:
        config["assembly_threads"]
    resources:
        mem = config["assembly_memory"],
        time= config["runtime"]["assembly"]
    shell:
        "spades.py --plasmid"
        " --threads {threads} "
        " --memory {resources.mem} "
        " -o {params.p[outdir]} "
        " -k {params.k}"
        " {params.p[preset]} "
        " {params.p[extra]} "
        " {params.p[inputs]} "
        " {params.p[longreads]} "
        " {params.p[skip_error_correction]} "
        " > {log} 2>&1 "
