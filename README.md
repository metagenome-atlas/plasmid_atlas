# Plasmid Atlas

This workflow is designed to run mataplasmidspades on the output of metagenome atlas

## Authors

* Silas Kieser (@silask)

## Dependencies

conda or vibrant

## Usage

This extension should be used after the running metagenome-atlas at least until the assembly step.


    ./plasmid_atlas atlas_working_directory path/to/store/databases {other snakemake arguments}

or directly call snakemake
    snakemake -d atlas_working_directory --config database_dir="path/to/store/databases"


For more detailed configuration see the `config/template_config.yaml`

# Cluster execution

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
