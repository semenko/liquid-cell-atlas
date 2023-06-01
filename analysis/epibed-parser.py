import pandas as pd
import pysam

def parse_methyl_string(methyl_string):
    state = methyl_string[0]
    count = 0
    result = []
    for char in methyl_string[1:]:
        if char.isdigit():
            count = count * 10 + int(char)
        else:
            result.append((state, count))
            state = char
            count = 0
    result.append((state, count))
    return result

def parse_epibed(file_path):
    bamfile = pysam.TabixFile(file_path)
    lines = bamfile.fetch()
    
    data = []
    for line in lines:
        fields = line.strip().split('\t')
        chrom, start, end, name, score, strand, methyl_string = fields[:7]
        data.append({
            'chrom': chrom,
            'start': int(start),
            'end': int(end),
            'name': name,
            'score': int(score),
            'strand': strand,
            'methylation': parse_methyl_string(methyl_string)
        })
    
    return pd.DataFrame(data)

df = parse_epibed('path/to/file.epibed.bgz')
