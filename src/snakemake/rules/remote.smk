from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule remote_http:
    input:
        #HTTP.remote("www.dropbox.com/s/dae6otw5honour0/Sequencing_summary.xlsx")
        HTTP.remote("dl.dropboxusercontent.com/s/dae6otw5honour0/Sequencing_summary.xlsx")
    output:
        #"out/remote/http/www.dropbox.com/s/dae6otw5honour0/Sequencing_summary.xlsx"
        "out/remote/http/dl.dropboxusercontent.com/s/dae6otw5honour0/Sequencing_summary.xlsx"
    shell:
        "mv {input} {output}"
