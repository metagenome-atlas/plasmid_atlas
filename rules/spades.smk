include: "sample_table.smk"

assembly_params=dict()
assembly_params['spades'] = {'meta':'--meta','normal':''}
JAVA_MEM_FRACTION = 0.85
MERGING_FLAGS = "ecct iterations=1"
MERGING_EXTEND2 = 50
MERGING_K = 62

ASSEMBLY_FRACTIONS=['R1','R2','me','se']
assembly_preprocessing_steps='QC.errorcorr.merged'


#SAMPLES= glob_wildcards(f"{{sample}}/plasmids/reads/{assembly_preprocessing_steps}_{ASSEMBLY_FRACTIONS[0]}.fastq.gz").sample


localrules: init_pre_assembly_processing
rule init_pre_assembly_processing:
    input:
        get_quality_controlled_reads
    output:
         temp(expand("{{sample}}/plasmids/reads/QC_{fraction}.fastq.gz",fraction= MULTIFILE_FRACTIONS))
    log:
        "{sample}/logs/plasmids/init.log"
    threads:
        1
    run:
        #make symlink
        assert len(input) == len(output), "Input and ouput files have not same number, can not create symlinks for all."
        for i in range(len(input)):
            os.symlink(os.path.abspath(input[i]),output[i])


rule error_correction:
    input:
        expand("{{sample}}/plasmids/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS)
    output:
        temp(expand("{{sample}}/plasmids/reads/{{previous_steps}}.errorcorr_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS))
    benchmark:
        "logs/benchmarks/plasmids/pre_process/{sample}_error_correction_{previous_steps}.txt"
    log:
        "{sample}/logs/plasmids/pre_process/error_correction_{previous_steps}.log"
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    params:
        inputs = lambda wc, input : io_params_for_tadpole(input),
        outputs = lambda wc, output: io_params_for_tadpole(output,key='out')
    threads:
        config.get("threads", 1)
    shell:
        """
        tadpole.sh -Xmx{resources.java_mem}G \
            prealloc=1 \
            {params.inputs} \
            {params.outputs} \
            mode=correct \
            threads={threads} \
            ecc=t ecco=t 2>> {log}
        """


rule merge_pairs:
    input:
        expand("{{sample}}/plasmids/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=['R1','R2'])
    output:
        temp(expand("{{sample}}/plasmids/reads/{{previous_steps}}.merged_{fraction}.fastq.gz",
            fraction=['R1','R2','me']))
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    conda:
        "../envs/bbmap.yaml"
    log:
        "{sample}/logs/plasmids/pre_process/merge_pairs_{previous_steps}.log"
    benchmark:
        "logs/benchmarks/plasmids/pre_process/merge_pairs_{previous_steps}/{sample}.txt"
    shadow:
        "shallow"
    params:
        kmer = config.get("merging_k", MERGING_K),
        extend2 = config.get("merging_extend2", MERGING_EXTEND2),
        flags = config.get("merging_flags", MERGING_FLAGS)
    shell:
        """
        bbmerge.sh -Xmx{resources.java_mem}G threads={threads} \
            in1={input[0]} in2={input[1]} \
            outmerged={output[2]} \
            outu={output[0]} outu2={output[1]} \
            {params.flags} k={params.kmer} \
            extend2={params.extend2} 2> {log}
        """
localrules: passtrough_se_merged
rule passtrough_se_merged:
    input:
        "{sample}/plasmids/reads/{previous_steps}_se.fastq.gz"
    output:
        temp("{sample}/plasmids/reads/{previous_steps}.merged_se.fastq.gz")
    shell:
        "cp {input} {output}"


#plasmid spades

def spades_parameters(wc,input):
    if not os.path.exists("{sample}/plasmids/spades/params.txt".format(sample=wc.sample)):

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


        params['preset'] =  ' --plasmid '+assembly_params['spades'][config['spades_preset']]
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
        expand("{{sample}}/plasmids/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
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
        "spades.py "
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


rule rename_plasmids:
    # standardizes header labels within contig FASTAs
    input:
        "{sample}/plasmids/spades/contigs.fasta"
    output:
        "{sample}/plasmids/plasmids.fasta.gz"
    conda:
        "../envs/bbmap.yaml"
    params:
        prefix="{sample}_circular"
    shell:
        "rename.sh in={input} out={output} ow=t prefix={params.prefix}"
