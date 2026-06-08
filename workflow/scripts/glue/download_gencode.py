if "snakemake" in locals():
    url = snakemake.params["url"]
    out_path = snakemake.output["gtf"]
else:
    url = (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
    )
    out_path = "data/reference/gencode.v38.primary_assembly.annotation.gtf.gz"

import gzip
import os
import urllib.request


tmp_path = f"{out_path}.download"

try:
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, tmp_path)

    if os.path.getsize(tmp_path) == 0:
        raise RuntimeError("Downloaded GENCODE annotation is empty.")

    with gzip.open(tmp_path, "rb") as handle:
        if not handle.read(1):
            raise RuntimeError("Downloaded GENCODE annotation has no content.")

    os.replace(tmp_path, out_path)
    print(f"Wrote GENCODE annotation to {out_path}")
finally:
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
