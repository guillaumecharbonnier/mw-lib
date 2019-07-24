rule conda_env_create_yaml:
    """
    Created:
        2018-08-17 00:54:36
    Aim:
        General purpose rule to install conda environment from yaml specs.
        Is mainly used for dev and debug.
        The specs are stored in snakemake tree so they can also be used by it with "--use-conda".
        Explicit listing of package in environment should help reproducing identical environments:
        https://conda.io/docs/user-guide/tasks/manage-environments.html
    Test:
        out/conda/env_create_yaml_list_explicit/r_component_analysis.txt
        out/conda/env_create_yaml_list_explicit/texlive.txt
        out/conda/env_create_yaml_list_explicit/snakemake.txt
    """
    input:
        #conda="opt/miniconda/bin/conda",
        yaml="src/snakemake/envs/{env}.yaml"
    output:
        #conda="opt/miniconda/envs/{env}/bin/conda",
        explicit="out/conda/env_create_yaml_list_explicit/{env}.txt"
    benchmark:
                 "out/conda/env_create_yaml_list_explicit/{env}.benchmark.tsv"
    shell:
        """
        conda env create --force --file {input.yaml} --name {wildcards.env} --verbose
        conda list --explicit --name {wildcards.env} > {output.explicit}
        cp {output.explicit} {output.explicit}.`date +%s`
        """


