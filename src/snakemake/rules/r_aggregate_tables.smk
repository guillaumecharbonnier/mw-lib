rule r_aggregate_tables:
    input:
        lambda wildcards: eval(mwconf['ids'][wildcards.table_list_id])
    output:
        "out/r/aggregate_tables/{table_list_id}"
    conda:
        "../envs/r.yaml"
    script:
        "../../r/scripts/aggregate_tables.R"

