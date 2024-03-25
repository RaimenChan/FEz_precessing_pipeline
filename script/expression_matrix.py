import os 
import pandas as pd
import argparse

##################################### arguments #####################################################

parser = argparse.ArgumentParser(description="This is a in-home script developed by Raiman (ruiwenchen@um.edu.mo) for FungiExpresZ RNAseq data processing pipeline.\n\
",\
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-alignmentstat_path", help="Path to the folder where sample_alignmentstat files will be stored.\n\n"
                                        )

parser.add_argument("-stringtie_path", help="Path to the folder where stringtie output sample.xls files will be stored.\n\n"
                                        )

parser.add_argument("-mapping_rate_t", type=float, help="Mapping rate threshold. \
                    Samples with both a mapping rate higher than the threshold and a mapping reads count higher than the reads threshold will be merged into the expression matrix.\
                    default = 70\n\n", default= 70
                                        )

parser.add_argument("-mapped_reads_t", type=int, help="Mapped reads threshold. \
                    Samples with both a mapping rate higher than the threshold and a mapping reads count higher than the reads threshold will be merged into the expression matrix.\
                    default = 0\n\n", default=0
                                        )

parser.add_argument("-o", help = "output path, default=./output", default='./output')


# read parameter
args = parser.parse_args()
alignmentstat_path = args.alignmentstat_path
stringtie_path = args.stringtie_path
mapping_rate_t = args.mapping_rate_t
mapped_reads_t = args.mapped_reads_t
output_path = args.o

# make output_folder
try:
    os.makedirs(output_path)
except:
    pass

## alignment stat
alignmentstats = os.listdir(alignmentstat_path)

mapped_readss = []
mapping_rates = []
names = []
for alignmentstat in alignmentstats:
    #alignmentstat_path of each file
    if(alignmentstat[-14:]=="_alignmentstat"):
        a_path2 = os.path.join(alignmentstat_path, alignmentstat)
        name = alignmentstat.split("_")[0]
        with open (a_path2, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            line6 = lines[6]
            mapped_reads = line6.split(" ")[0]
            mapping_rate = line6.split("(")[1].split(" :")[0]
            mapping_rate = float(mapping_rate.strip('%'))
            
        names.append(name)   
        mapped_readss.append(mapped_reads)
        mapping_rates.append(mapping_rate)
            
alignment_df = pd.DataFrame()
alignment_df['sample'] = names
alignment_df['mapped_reads'] = mapped_readss
alignment_df['mapped_reads'] = alignment_df['mapped_reads'].astype(int)
alignment_df['mapping_rate'] = mapping_rates
alignment_df['mapping_rate'] = alignment_df['mapping_rate'].astype(float)
alignment_df = alignment_df.sort_values(['mapping_rate', 'mapped_reads'], ascending=[False, False])

alignment_op = os.path.join(output_path, "alignment.xlsx")
alignment_df.to_excel(alignment_op, index=False)


# cutoff_samples
cf_samples = alignment_df[ (alignment_df['mapping_rate']>mapping_rate_t) & (alignment_df['mapped_reads']>mapped_reads_t)]['sample'].tolist()


## merge good sample (cf_samples) to expression matrix
flag = 0 
for i,sample in enumerate(cf_samples):
    #sample_path of each file
    s_path2 = os.path.join(stringtie_path, "".join([sample, ".xls"]))
    try:
        sample_df = pd.read_csv(s_path2, sep='\t', comment="#", header=None)
        sample_df = sample_df[sample_df[2]=='transcript']
        sample_df['ID'] = sample_df[8].apply(lambda x:x.split("\";")[0].split("\"")[1])
        sample_df[sample] = sample_df[8].apply(lambda x:x.split("\";")[3].split("\"")[1])
        sample_df = sample_df[['ID', sample]]
        if(flag==0):
            print("i=%d, mergeing %s to expression matrix"%(i, sample))
            samples_df = sample_df
            flag=1
        else:
            print("i=%d, mergeing %s to expression matrix"%(i, sample))
            samples_df = pd.merge(samples_df, sample_df, on='ID', how='outer')
    except:
        print("i=%d, could not open %s. skipping this sample."%(i, s_path2))


samples_df = samples_df.fillna(0)
expression_op = os.path.join(output_path, "expression_matrix.tsv")
samples_df.to_csv(expression_op, sep='\t', index=None)


print("end")



