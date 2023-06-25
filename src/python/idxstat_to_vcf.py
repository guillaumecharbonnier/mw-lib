#!/usr/bin/env python


#GC I have a file produced from samtools idxstats looking like this:
# Contig1_1:959510_REF_A  41      0       0
# Contig1_1:959510_ALT_AAGAC      45      2       0
# Contig2_1:1013466_REF_T 41      4       0
# Contig2_1:1013466_ALT_TA        42      2       0
# Contig3_1:2416398_REF_G 41      0       0
# Contig3_1:2416398_ALT_GGGACA    46      1       0
# Contig4_1:3625542_REF_T 41      0       0
# Contig4_1:3625542_ALT_TG        42      2       0
# 
# I would like to write a script that convert it into a vcf file where the CHROM is taken from the number after the first underscore, the POS from the number after the first :, ID from the Contig before the first underscore, REF from the nucleotidic sequence on rows containing "REF", ALT from the nucleotidic sequence on rows containing "ALT".  
# For the integer columns of the idxstat file, I only need to extract the second one, as NREF for the REF row and as NALT for the ALT row.


# import re

# with open('idxstats.txt', 'r') as f:
#     lines = f.readlines()

# with open('output.vcf', 'w') as f:
#     for line in lines:
#         if line.startswith('*'):
#             continue
#         fields = line.strip().split()
#         contig, chrom, pos, ref_alt = re.match(r'(\w+)_(\d+):(\d+)_(\w+)', fields[0]).groups()
#         nref_nalt = fields[2]
#         if 'REF' in ref_alt:
#             ref = ref_alt.split('_')[-1]
#             f.write(f'{chrom}\t{pos}\t{contig}\t{ref}\t.\t.\t.\tNREF={nref_nalt}')
#         else:
#             alt = ref_alt.split('_')[-1]
#             f.write(f'\t{alt}\t.\t.\tNALT={nref_nalt}\n')


#GC I still get the error. I believe i may have some other rows that do not match the regex. Could you give me code to highlight such rows ? Currently the error does not higlight them

# Sure! Here is a modified version of the script that will print out any lines that do not match the regular expression:

# import re

# with open('idxstats.txt', 'r') as f:
#     lines = f.readlines()

# with open('output.vcf', 'w') as f:
#     for line in lines:
#         if line.startswith('*'):
#             continue
#         fields = line.strip().split()
#         match = re.match(r'(\w+)_(\d+):(\d+)_(\w+)', fields[0])
#         if not match:
#             print(f'Line does not match regular expression: {line}')
#             continue
#         contig, chrom, pos, ref_alt = match.groups()
#         nref_nalt = fields[2]
#         if 'REF' in ref_alt:
#             ref = ref_alt.split('_')[-1]
#             f.write(f'{chrom}\t{pos}\t{contig}\t{ref}\t.\t.\t.\tNREF={nref_nalt}')
#         else:
#             alt = ref_alt.split('_')[-1]
#             f.write(f'\t{alt}\t.\t.\tNALT={nref_nalt}\n')

# This script will print out any lines from your idxstats file that do not match the regular expression `r'(\w+)_(\d+):(\d+)_(\w+)'`. You can then use this information to identify any issues with the format of your idxstats file and adjust the script accordingly.

#GC Great, I get the error now. It is because the CHROM can be either an integer or something more complex. I think it would be better to capture the whole string between the first _ and : as chrom because there are also special chrom that may contains _. 

# Sure! Here is an updated version of the script that captures the entire string between the first underscore (`_`) and the first colon (`:`) as the `chrom` field:

# import re

# with open('idxstats.txt', 'r') as f:
#     lines = f.readlines()

# with open('output.vcf', 'w') as f:
#     for line in lines:
#         if line.startswith('*'):
#             continue
#         fields = line.strip().split()
#         match = re.match(r'(\w+)_(.+):(\d+)_(\w+)', fields[0])
#         if not match:
#             print(f'Line does not match regular expression: {line}')
#             continue
#         contig, chrom, pos, ref_alt = match.groups()
#         nref_nalt = fields[2]
#         if 'REF' in ref_alt:
#             ref = ref_alt.split('_')[-1]
#             f.write(f'{chrom}\t{pos}\t{contig}\t{ref}\t.\t.\t.\tNREF={nref_nalt}')
#         else:
#             alt = ref_alt.split('_')[-1]
#             f.write(f'\t{alt}\t.\t.\tNALT={nref_nalt}\n')

# This updated script should now correctly handle `chrom` values that contain underscores (`_`) or other special characters. Let me know if you have any further issues or questions.

#GC Great, I can now see an output file. I have changed my mind and I would like to replace NALT by NV and NREF should not be used. Instead, I need NR which should be NALT+NREF

# Sure! Here is an updated version of the script that replaces `NALT` with `NV` and calculates `NR` as the sum of `NALT` and `NREF`:

# import re

# with open('idxstats.txt', 'r') as f:
#     lines = f.readlines()

# with open('output.vcf', 'w') as f:
#     nref = 0
#     for line in lines:
#         if line.startswith('*'):
#             continue
#         fields = line.strip().split()
#         match = re.match(r'(\w+)_(.+):(\d+)_(\w+)', fields[0])
#         if not match:
#             print(f'Line does not match regular expression: {line}')
#             continue
#         contig, chrom, pos, ref_alt = match.groups()
#         nref_nalt = int(fields[2])
#         if 'REF' in ref_alt:
#             ref = ref_alt.split('_')[-1]
#             nref = nref_nalt
#             f.write(f'{chrom}\t{pos}\t{contig}\t{ref}\t.\t.\t.\t')
#         else:
#             alt = ref_alt.split('_')[-1]
#             nr = nref + nref_nalt
#             f.write(f'{alt}\t.\t.\tNV={nref_nalt};NR={nr}\n')

# This updated script should now produce a VCF file with the desired format. Let me know if you have any further issues or questions.

#GC Nice! Here is a few rows of the output:
# 1       959510  Contig1 A       .       .       .       AAGAC   .       .       NV=2;NR=2
# 1       1013466 Contig2 T       .       .       .       TA      .       .       NV=2;NR=6
# 1       2416398 Contig3 G       .       .       .       GGGACA  .       .       NV=1;NR=1
# 
# What are those dots columns? I would like to match the format of a VCF file:
# 1	CHROM
# 2	POS
# 3	ID
# 4	REF
# 5	ALT
# 6	QUAL
# 7	FILTER
# 8	INFO
# 
# For QUAL and FILTER we can put a dot. 

# The dots in the output file represent missing values for the `QUAL` and `FILTER` columns of the VCF file. In a VCF file, a dot (`.`) is used to indicate that a value is missing or not applicable.

# Here is an updated version of the script that produces a VCF file with the desired column headers and format:

import re

with open('idxstats.txt', 'r') as f:
    lines = f.readlines()

with open('output.vcf', 'w') as f:
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    nref = 0
    for line in lines:
        if line.startswith('*'):
            continue
        fields = line.strip().split()
        match = re.match(r'(\w+)_(.+):(\d+)_(\w+)', fields[0])
        if not match:
            print(f'Line does not match regular expression: {line}')
            continue
        contig, chrom, pos, ref_alt = match.groups()
        nref_nalt = int(fields[2])
        if 'REF' in ref_alt:
            ref = ref_alt.split('_')[-1]
            nref = nref_nalt
            f.write(f'{chrom}\t{pos}\t{contig}\t{ref}\t')
        else:
            alt = ref_alt.split('_')[-1]
            nr = nref + nref_nalt
            f.write(f'{alt}\t.\t.\tNV={nref_nalt};NR={nr}\n')

# This updated script should now produce a VCF file with the desired column headers and format. Let me know if you have any further issues or questions.