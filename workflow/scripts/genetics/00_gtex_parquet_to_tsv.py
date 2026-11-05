# (C) Sergio Alejandro Gómez
# Aim: Convert the GTEx v10 kidney-cortex significant-pairs parquet to a gzip TSV that R can read (R has no parquet
#   reader here). Lossless: row order and all columns kept, floats written with 17 significant digits, gzip mtime = 0.
# Inputs: Kidney_Cortex.v10.eQTLs.signif_pairs.parquet
# Outputs: gzip TSV of the same table
# Usage: python workflow/scripts/genetics/00_gtex_parquet_to_tsv.py --parquet=data/genetics/Kidney_Cortex.v10.eQTLs.signif_pairs.parquet --out=results/genetics/qtl_resources/gtex_kidney_signif_pairs.tsv.gz
import sys, io, gzip
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import pyarrow.parquet as pq

args = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
table = pq.read_table(args['parquet']).to_pandas()
with open(args['out'], 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, compresslevel=6, mtime=0) as gz, \
        io.TextIOWrapper(gz, encoding='utf-8', newline='\n') as fh:
    table.to_csv(fh, sep='\t', index=False, float_format='%.17g', lineterminator='\n')
print(f"wrote {args['out']} ({len(table)} rows, {table.shape[1]} columns)")
