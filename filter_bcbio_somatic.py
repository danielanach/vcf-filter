from pysam import VariantFile
import pandas as pd
import string

def filter_bcbio_somatic(in_vcf_name,
                         out_vcf_name,
                         max_maf=0.001,
                         max_af=0.9,
                         loh_df=None):

    vcf_in = VariantFile(in_vcf_name)

    out_vcf = open(out_vcf_name,mode='w')
    out_vcf.write(vcf_in.header.__str__())

    df = pd.read_csv('/Users/d/work/dcis_lcm_exome/data/cancer_genes.tsv',sep='\t')
    cancer_gene_list = [g.strip() for g in df['Gene Symbol']]

    for rec in vcf_in:

        write_var = True
        if 'max_aaf_all' in rec.info.keys():
            # Found in population databases
            if is_frequent(rec,max_maf):

                if not in_cancer_gene(rec,cancer_gene_list):
                    continue

                if not is_in_clinvar(rec) and not is_cosmic_variant(rec):
                    continue

                if is_in_clinvar(rec) and is_benign_in_clinvar(rec):
                    continue

                if is_clonal_and_not_loh(rec,max_af,loh_df):
                    continue
                else:
                    out_vcf.write(rec.__str__())
        else:
            if is_clonal_and_not_loh(rec,max_af,loh_df):
                continue
            else:
                try:
                    out_vcf.write(rec.__str__())
                except:
                    s = rec.info['dbNSFP_clinvar_trait']
                    printable = set(string.printable)
                    filter(lambda x: x in printable, s)
                    rec.info['dbNSFP_clinvar_trait'] = ''.join(filter(lambda x: x in printable, s))
                    out_vcf.write(rec.__str__())

def in_cancer_gene(rec,cancer_gene_list):
    found_t = False
    for a in rec.info.get('ANN'):
        if a.split('|')[5].strip() == 'transcript':
            found_t = True
            ann = a
    if not found_t:
        ann = rec.info.get('ANN')[0]
    ann = ann.split('|')
    gene = ann[3]

    return gene in cancer_gene_list

def is_clonal_and_not_loh(rec,max_af,loh_df):

    likely_germ = False

    if (rec.info['AF'][0] >= max_af) and (loh_df is not None):
            if not var_in_loh(loh_df,rec.contig,rec.start,rec.stop):
                likely_germ = True

    return likely_germ

def var_in_loh(loh_df,chrom,start,stop):

    in_loh = False
    for i in loh_df.index:
        reg = loh_df.loc[i]
        if reg['chromosome'] == chrom:
            if start >= reg['start'] and stop <= reg['end']:
                break
    return in_loh

def is_frequent(rec,max_maf):

    frequent = False

    if 'max_aaf_all' in rec.info.keys():
        if isinstance(rec.info['max_aaf_all'], float):
            max_aaf_all = float(rec.info['max_aaf_all'])
        else:
            max_aaf_all = float(rec.info['max_aaf_all'][0])

        if max_aaf_all > max_maf:
            frequent = True

    return frequent

def is_in_clinvar(rec):
    if ('clinvar_sig' in rec.info.keys()) or ('dbNSFP_clinvar_clnsig' in rec.info.keys()):
        return True
    else:
        return False

def is_benign_in_clinvar(rec):

    benign = False
    if 'clinvar_sig' in rec.info.keys():
        sig = ''.join(rec.info['clinvar_sig']).lower()
        if 'benign' in sig:
            benign = True
    if 'dbNSFP_clinvar_clnsig' in rec.info.keys():
        if ('2' in ''.join(rec.info['dbNSFP_clinvar_clnsig'])) or ('3' in ''.join(rec.info['dbNSFP_clinvar_clnsig'])):
            benign = True

    return benign

def is_cosmic_variant(rec):

    cosmic = False

    if 'cosmic_ids' in rec.info.keys():
        if not 'COSN' in "".join(rec.info['cosmic_ids']):
            cosmic = True

    return cosmic
