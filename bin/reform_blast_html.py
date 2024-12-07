#!/usr/bin/env python3
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description='Reform BLAST Output HTML Page.')
parser.add_argument('--tsv', type=str, dest='tsv_doc', required=True,
                    help='BLAST output TSV file')
parser.add_argument('--html', type=str, dest='html_in', required=True,
                    help='BLAST output HTML file')

def get_accession_dict(tsv_doc):
    try:
        blast_df = pd.read_csv(tsv_doc, comment="#", sep="\t", header=None)
        if not blast_df.empty:
            blast_df.columns =['acc_ver', 'identity', 'bit_score', 'length', 'evalue', 'title']
            return dict(zip(blast_df['acc_ver'], blast_df['identity']))
        else:
            return None
    except:
        return None

def main():
    args = parser.parse_args()
    html_in = args.html_in
    html_out = html_in.split(".")[0] + "_identity.html"
    acc_ver_dict = get_accession_dict(args.tsv_doc)
    if acc_ver_dict != None:
        with open(html_in) as f, open(html_out,'w') as of:
            for line in f:
                if line.startswith("<b>BLASTN"):
                    of.write(line.replace("BLASTN","BLAST ®"))
                elif "<b>Query=</b>" in line:
                    f.readline()
                    of.write(f.readline().replace("E","Identity"))
                    of.write(f.readline().replace("Value","  (Percent%)"))
                elif line.startswith('<a title='):
                    items = line.strip().split("</a>")
                    items.pop()
                    newline = '</a>'.join(items)
                    m = re.match("<a title=\"Show report for (.*[0-9]+\\.[0-9]+)\".*",items[0])
                    of.write("%s\n"%(newline + "</a>\t" + str(acc_ver_dict[m.group(1)])))
                else:
                    of.write(line)
    else:
        with open(html_in) as f, open(html_out,'w') as of:
            for line in f:
                of.write(line)

if __name__ == "__main__":
    main()