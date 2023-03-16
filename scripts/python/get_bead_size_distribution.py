import pandas as pd
import numpy as np
import argparse
import glob
import tqdm
from collections import defaultdict
from matplotlib import cm
import matplotlib
matplotlib.use('Agg')

'''
Profiles the proportion of clusters and proportion of reads within various cluster size categories. Considers DPM and BPM reads seperately (based on the input parameter readtype).
'''


def main():
    args = parse_arguments()
    pattern = args.directory + "/*" + args.pattern
    files = glob.glob(pattern)
    cluster_counts = []
    read_counts = []
    for f in files:
       df1, df2 = count_statistics(f, args.readtype) 
       cluster_counts.append(df1)
       read_counts.append(df2)
    cluster_df = pd.concat(cluster_counts, axis=1).transpose()
    read_df = pd.concat(read_counts, axis=1).transpose()
    cluster_fig = plot_profile(cluster_df)
    cluster_fig.savefig(args.directory + '/'+ args.readtype + '_cluster_distribution.pdf', bbox_inches='tight')
    read_fig = plot_profile(read_df)
    read_fig.savefig(args.directory + '/' + args.readtype + '_read_distribution.pdf', bbox_inches='tight')

def count_statistics(filename, readtype):
    '''
    Loop over all clusters in a clusterfile and bin read/cluster counts based on the size of the cluster.

    Note:
        size of cluster only includes the number of reads corresponding to the read type being counted

    Args:
        filename(str): path to clusterfile)
        readtype(str): DPM or BPM
    '''
    bins = np.array([0,1,5,10,20, 50,100,200]) #bins with data are 1-8
    cluster_counts = defaultdict(int)
    read_counts = defaultdict(int)
    for bin in np.arange(1,len(bins)+1): #initialize all bins in case any end up being empty categories
        cluster_counts[bin] = 0
        read_counts[bin] = 0
    with open(filename, 'r') as f:
        for line in tqdm.tqdm(f):
            if readtype in line:
                reads = line.split('\t',1)[1] #ignore barcode
                counts = reads.count(readtype)
                bin = np.digitize(counts, bins, right=True)
                cluster_counts[bin] +=1
                read_counts[bin] += counts
    df_cluster_counts = pd.DataFrame.from_dict(cluster_counts, orient='index')
    df_cluster_counts.columns = [filename.rsplit('/',1)[-1].rsplit('.clusters')[0]]
    df_read_counts = pd.DataFrame.from_dict(read_counts, orient='index')
    df_read_counts.columns = [filename.rsplit('/',1)[-1].rsplit('.clusters')[0]]
    return df_cluster_counts, df_read_counts

def plot_profile(df):
    '''
    Plot the proportion of reads within each size category as a stacked bar graph
    Args:
        df(dataframe): binned counts 
    '''
    columns = ['1', '2-5', '6-10', '11-20', '21-50', '51-100', '101-200', '201+']
    df.columns = columns
    df = df.div(df.sum(axis=1), axis=0) 
    plot = df.plot(kind='bar',stacked=True, ylabel='Proportion', cmap = cm.get_cmap('Dark2'))
    plot.legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
    plot.set_ylabel = 'Proportion'
    return plot.get_figure()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generate the cluster size distribution plot.')

    parser.add_argument('--directory',
                        metavar = "FILE",
                        action = "store",
                        help = "The directory with clusterfiles")
    parser.add_argument('--pattern',
                        metavar = "FILE",
                        action = "store",
                        help = "The pattern for file names")
    parser.add_argument('--readtype',
                        action="store",
                        help = "The read type to count")
    return parser.parse_args()

if __name__ == "__main__":
    main()
            












