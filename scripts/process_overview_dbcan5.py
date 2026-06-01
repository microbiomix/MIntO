#!/usr/bin/env python

import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="Parse run_dbCANv5 overview.tsv output"
)
parser.add_argument('-i', '--input', help='Input overview.tsv')
parser.add_argument('-o', '--output', help='Output tsv path prefix, - for stdout')
parser.add_argument('--format', choices=['minto', 'long'], default='minto',
                    help='Output format comma-separated compact (minto) or long')

args = parser.parse_args()

def make_unique_words(x):
    return(','.join(set(x.split(","))))
    
def make_unique_domains(x):
    return(','.join(sorted(list(set(x.split("|"))))))

input_file = args.input
output_prefix = args.output

colnames = ['ID', 'EC_nums', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', 'NumberofTools', 'RecommendResults', 'Substrate']
df = pd.read_csv(input_file, header=0, names=colnames, sep = "\t")

# process only those with min 2 tools agree
df = df[df.NumberofTools > 1]

# remove some cols
df.drop(columns = ['dbCAN_hmm','dbCAN_sub','DIAMOND', 'NumberofTools'], inplace=True)

if (args.format == "minto"):
    # replace hit separators with comma
    df['EC_nums'] = df['EC_nums'].str.replace(r':\d+', '', regex=True).str.replace(";", ",")
    df['Substrate'] = df['Substrate'].str.replace(";", ",")
    df['RecommendResults'] = df['RecommendResults'].transform({'RecommendResults': make_unique_domains})
    df['EC_nums'] = df['EC_nums'].transform({'EC_nums': make_unique_domains})
    df['Substrate'] = df['Substrate'].transform({'Substrate': make_unique_words})
    # extract the binding domains
    df['dbCAN.binding_module'] = df['RecommendResults'].str.replace(r',[CG][ETH].*', '', regex=True).str.extract("(CBM\d+.*),?G?")
    # remove the eCAMI sub-clustering numbers
    df['dbCAN.binding_module'] = df['dbCAN.binding_module'].str.replace(r'_e\d+', '', regex=True)
    df['dbCAN.subfamily'] = df['RecommendResults'].str.replace(r'_e\d+', '', regex=True)
    # extract the cazymes
    df['dbCAN.subfamily'] = df['dbCAN.subfamily'].str.replace(r'CBM\d+,?', '', regex=True).str.replace(r',$', '', regex=True)
    # fill NAN and rename
    df = df.fillna("-").rename(columns={"EC_nums": "dbCAN.EC", "RecommendResults": "dbCAN.RecommendResults", "Substrate": "dbCAN.Substrate"})
else:
    # separate the multi-domain hits into list-in-cells
    df["dbCAN.hit"]=df["RecommendResults"].str.split("|")
    df["dbCAN.EC"]=df["EC_nums"].str.split("|")
    # separate list-in-cell into rows for dbCAN.hit
    df1 = df.explode("dbCAN.hit").reset_index(drop=True)
    df1["rownum"] = df1.groupby("ID")["dbCAN.hit"].cumcount().add(1)
    df1.drop(columns = ['EC_nums','RecommendResults','dbCAN.EC'], inplace=True)
    # separate list-in-cell into rows for EC numbers
    df2 = df.explode("dbCAN.EC").reset_index(drop=True)
    df2["rownum"] = df2.groupby("ID")["dbCAN.EC"].cumcount().add(1)
    df2 = df2[["ID", "dbCAN.EC", "rownum"]]
    # in theory each domain separated by | have an EC number so df1 and df2 should have corresponding amount of rows
    df = df1.merge(df2, on = ['ID', 'rownum'], how = "left")
    df = df.fillna("-")
    df.drop(columns = ['rownum'], inplace=True)
    # remove duplicates
    df = df.drop_duplicates()
    # CBM will be bindingmodule, rest cazyme; there can be multiple of both
    df["dbCAN.binding_module"] = '-'
    df.loc[df['dbCAN.hit'].str.contains("CBM"), 'dbCAN.binding_module'] = df['dbCAN.hit']
    df["dbCAN.cazyme"] = '-'
    df.loc[df['dbCAN.hit'].str.contains("CBM") == False, 'dbCAN.cazyme'] = df['dbCAN.hit']
    # remove the eCAMI sub-clustering numbers, also from CBM
    df['dbCAN.subfamily'] = df['dbCAN.cazyme'].str.replace(r'_e\d+', '', regex=True)
    df['dbCAN.binding_module'] = df['dbCAN.binding_module'].str.replace(r'_e\d+', '', regex=True)
    # strip EC
    df['dbCAN.EC'] = df['dbCAN.EC'].str.replace(r':\d+', '', regex=True).str.replace(";", ",")
    df['Substrate'] = df['Substrate'].str.replace(";", ",")
    # rename cols
    df = df.rename(columns={"dbCAN.hit": "dbCAN.RecommendResults", "Substrate": "dbCAN.Substrate"})

if (output_prefix != "-"):
    # output subfamily file and substrate separately
    output_file = f"{output_prefix}.domains.tsv"
    df.to_csv(output_file, columns = ["ID", "dbCAN.RecommendResults", "dbCAN.binding_module", "dbCAN.subfamily", "dbCAN.EC"], sep = "\t", index=False)

    output_substr = f"{output_prefix}.substrate.tsv"
    df = df[['ID', 'dbCAN.Substrate']]
    df[df['dbCAN.Substrate'] != "-"].drop_duplicates().to_csv(output_substr, sep = "\t", index=False)
else:
    # stdout and full df
    df.to_csv(sys.stdout, sep = "\t", index=False)
