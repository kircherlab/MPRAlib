import pandas as pd
import os


def chromosome_map() -> pd.DataFrame:
    base_path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_path, "data", "hg19.chromAlias.txt")
    df = pd.read_csv(file_path, sep="\t", header=None, comment="#", dtype="category")
    file_path = os.path.join(base_path, "data", "hg38.chromAlias.txt")
    df = pd.concat(
        [
            df,
            pd.read_csv(
                file_path, sep="\t", header=None, comment="#", dtype="category"
            ),
        ],
        ignore_index=True,
    )
    df.columns = ["ucsc", "assembly", "genbank", "refseq"]
    return df
