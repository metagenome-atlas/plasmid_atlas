
import os,sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
from workflow_scripts import *
from snakemake.utils import validate


configfile: os.path.join(os.path.dirname(workflow.snakefile),'config/template_config.yaml'), 'config.yaml'
validate(config, "config/config.schema.yaml")


include: "rules/spades.smk"


rule all:
    input:
        expand("{sample}/plasmids/spades/contigs.fasta",sample=get_all_samples(config["sampletable"])[:2])



for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]
    if not "time" in r.resources:
        r.resources["time"]=config["time"]


#
