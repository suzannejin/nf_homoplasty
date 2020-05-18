#!/nfs/software/as/el7.2/EasyBuild/CRG/software/Python/3.7.2-GCCcore-8.2.0/bin/python3.7

# Python version 3.7

__description__ = '''

Normalize raw data.

'''


if __name__ == '__main__':

    import argparse
    import pandas as pd
    import sys

    app = argparse.ArgumentParser(description=__description__)
    app.add_argument('-i', type=str, help="Input data")
    app.add_argument('-o', type=str, help="Output data")
    app.add_argument('-norm', type=str, nargs='+', help="Columns to be normalized (divide by len & nseq")
    args = app.parse_args()

    # Read input
    df = pd.read_csv(args.i, header=0)

    # Normalize columns
    df[args.norm] = df[args.norm].div(df.nseq, axis = 0).div(df.len, axis = 0)

    # Write output
    df.to_csv(args.o, sep=",", na_rep='NA', header=True, index=False)
