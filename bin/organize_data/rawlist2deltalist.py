#!/nfs/software/as/el7.2/EasyBuild/CRG/software/Python/3.7.2-GCCcore-8.2.0/bin/python3.7

# Python version 3.7

__description__ = '''

Calculate the delta from the raw data.

'''


def rawlist2deltalist(df, maintain, deltaby, metrics, norm, decimals, log=None):

    import datetime
    import numpy as np
    import pandas as pd

    # Normalize if needed
    if norm:
        df[norm] = df[norm].div(df.nseq, axis = 0).div(df.len, axis = 0)

    # Initialize delta dataframe
    variables = maintain[:]
    variables.append(deltaby+"1")
    variables.append(deltaby+"2")
    header = variables + metrics  # Dataframe header
    delta = pd.DataFrame(columns = header)

    # Get delta group by group
    g = list(df.groupby(maintain)) 
    #print(len(g))
    for element in g:
        comb = element[1]   # Reduced dataframe grouby typs
        keys = comb[deltaby].tolist()
        nkeys = len(keys)
        for i in range(nkeys):

            # Get metrics and vectorize
            tmp = comb[metrics].to_numpy()
            values1 = tmp[i]
            values2 = np.delete(tmp, i, 0)
            # Calculate delta
            vals_delta = values1 - values2
            
            # Get fixed names
            variables = maintain[:]
            variables.append(deltaby)
            names1 = np.array( comb[variables].to_numpy()[[i]].tolist() * (nkeys-1) )
            names2 = comb[deltaby].to_numpy()
            names2 = np.delete(names2, i, 0).reshape(-1,1)
            names = np.append(names1, names2, axis=1)

            # Merge names + values
            #vals_delta = vals_delta.round(decimals)
            np_delta = np.append(names, vals_delta, axis=1)
            df_delta = pd.DataFrame(np_delta, columns = header)

            # Concatenate rows
            delta = pd.concat( [delta, df_delta] )
    
        # Print Process - log
        if log: 
            o = open(log, 'a')
            date = datetime.datetime.now()
            o.write( str(element[0]) + " " + str(date) + "\n" )

    if log:
        o.close()

    return(delta)


if __name__ == '__main__':

    import argparse
    import pandas as pd
    import sys

    app = argparse.ArgumentParser(description=__description__)
    app.add_argument('-raw', type=str, help="Raw data in list")
    app.add_argument('-delta', type=str, help="Output delta in list")
    app.add_argument('-maintain', type=str, nargs='+', help="Columns to be maintained while calculating delta")
    app.add_argument('-deltaby', type=str, help="Column considered to calculate the delta")
    app.add_argument('-metrics', type=str, nargs='+', help="The metrics used. Eg.: tc sp col homo whomo whomo2 len ngap ngap2 SSPA")
    app.add_argument('-norm', type=str, nargs='+', default=None, help="Metrics to be normalized (divide by nseq and len)")
    app.add_argument('-decimals', type=int, default=3)
    app.add_argument('-log', type=str, default=None, help="Log file")
    args = app.parse_args()

    # Print log
    if args.log:
        o = open(args.log, 'w')
        o.write("Starting converting raw list to delta list with the following command line:\n")
        o.write(" ".join(sys.argv) + "\n")
        o.close()

    # Convert raw list to delta list
    raw = pd.read_csv(args.raw, header=0)
    delta = rawlist2deltalist(raw, args.maintain, args.deltaby, args.metrics, args.norm, args.decimals, args.log)
    delta.to_csv(args.delta, sep=",", na_rep='NA', header=True, index=False)
