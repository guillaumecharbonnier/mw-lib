import sys
import mysql.connector
import pandas as pd

# Get the gene name and genome from the command line arguments
gene_name = sys.argv[1]
genome = sys.argv[2]
outfile = sys.argv[3]

# Connect to the UCSC MySQL server
cnx = mysql.connector.connect(
    user="genome", host="genome-mysql.soe.ucsc.edu", database=genome
)

# Create a cursor object
cursor = cnx.cursor()

# Execute the query
query = f"""
SELECT 
    kg.chrom,
    kg.txStart,
    kg.txEnd,
    kg.name,
    0 as score,
    kg.strand,
    kg.cdsStart,
    kg.cdsEnd,
    '0,0,0' as itemRgb,
    kg.exonCount,
    kg.exonStarts,
    kg.exonEnds
FROM 
    knownGene as kg
JOIN 
    kgXref as xref ON kg.name = xref.kgID
JOIN
    {genome}.knownCanonical as canonical ON kg.name = canonical.transcript
WHERE 
    xref.geneSymbol = '{gene_name}'
"""
cursor.execute(query)

# Fetch the results into a pandas DataFrame
df = pd.DataFrame(cursor.fetchall(), columns=[i[0] for i in cursor.description])

# Close the cursor and connection
cursor.close()
cnx.close()

# Combine the name and " gene=" + gene_name into the 4th column
# df["name"] = df["name"] + " gene=" + gene_name

# Process the exonStarts and exonEnds to compute exon sizes and relative starts
df["exonStarts"] = (
    df["exonStarts"].apply(lambda x: x.decode("utf-8")).str.rstrip(",").str.split(",")
)
df["exonEnds"] = (
    df["exonEnds"].apply(lambda x: x.decode("utf-8")).str.rstrip(",").str.split(",")
)
df["blockSizes"] = df.apply(
    lambda row: [
        int(end) - int(start) for start, end in zip(row["exonStarts"], row["exonEnds"])
    ],
    axis=1,
)
df["blockStarts"] = df.apply(
    lambda row: [int(start) - row["txStart"] for start in row["exonStarts"]], axis=1
)

# Convert lists to comma-separated strings
df["blockSizes"] = df["blockSizes"].apply(lambda x: ",".join(map(str, x)) + ",")
df["blockStarts"] = df["blockStarts"].apply(lambda x: ",".join(map(str, x)) + ",")

# Write the DataFrame to a BED file, excluding the 'exonStarts' and 'exonEnds' columns
df.drop(columns=["exonStarts", "exonEnds"]).to_csv(
    f"{outfile}", sep="\t", header=False, index=False
)
