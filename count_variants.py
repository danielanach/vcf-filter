import pandas as pd
from pysam import VariantFile

def count_variants(vcf_list, sample_list):
    """Counts variants in vcf file and outputs summary dataframe.

    Parameters
    ----------
    in_vcf : list(str)
      list of VCF file paths.

    out_vcf : list(str)
      list of sample names.

    Returns:

    pd DataFrame where columns are SNV and indel Counts and rows are samples.

    """

    df_lst = []

    sample_vcf_dct = dict(zip(sample_list,vcf_list))

    for s in sample_vcf_dct.keys():

        vcf_in = sample_vcf_dct[s]
        vcf = VariantFile(vcf_in)

        snv = 0
        indel = 0

        for rec in vcf:

            ref_len = len(rec.ref)

            for a in rec.alts:
                if len(a) > 1 or ref_len > 1:
                    indel +=1
                else:
                    snv +=1

        df_lst.append([s,snv,indel])

    out_df = pd.DataFrame(df_lst, columns=['sample','snvs','indels'])

    return out_df
