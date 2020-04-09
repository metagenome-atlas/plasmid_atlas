
import os,sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
from workflow_scripts import *
from snakemake.utils import validate


configfile: os.path.join(os.path.dirname(workflow.snakefile),'config/template_config.yaml')
validate(config, "config/config.schema.yaml")




include: "rules/spades.smk"
include: "rules/dereplicate.smk"


rule all:
    input:
        plasmids=expand("{sample}/plasmids/plasmids.fasta.gz",sample=SAMPLES),
        dereplicated="plasmids/plasmids.deduplicated.fasta.gz"




for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]
    if not "time" in r.resources:
        r.resources["time"]=config["time"]


#
