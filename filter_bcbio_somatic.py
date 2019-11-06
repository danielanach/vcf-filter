from pysam import VariantFile

def filter_bcbio_somatic(in_vcf_name,
                         out_vcf_name,
                         max_maf,
                         max_af=0.9,
                         loh_sites=[],
                         non_path_lst=['benign'],
                         loh_df=None):

    """Filters the bcbio output ensemble 'somatic' VCF which had somatic prioritization via gemini.

    https://bcbio.wordpress.com/2015/03/05/cancerval/

    Parameters
    ----------
    in_vcf : str
        Path to input VCF file.
    out_vcf : str
        Path to output VCF file.
    max_maf : float
        Max minor allelic fraction in any population to consider if not in COSMIC or pathogenic in clinvar.
    max_af : float
        Maximum variant allelic fraction to be considered somatic, only used when paired with loh_df
    non_path_lst : lst(str)
        List of clinvar labels to use as "non-pathogenic"

    """

    vcf_in = VariantFile(in_vcf_name)

    out_vcf = open(out_vcf_name,mode='w')
    out_vcf.write(vcf_in.header.__str__())

    for rec in vcf_in:

        info_dct = rec.info

        if (info_dct['AF'][0] >= max_af) and (loh_df is not None):
            if not var_in_loh(loh_df,rec.contig,rec.start,rec.stop):
                continue

        if 'max_aaf_all' in info_dct.keys():

            if isinstance(info_dct['max_aaf_all'], float):
                max_aaf_all = float(info_dct['max_aaf_all'])
            else:
                max_aaf_all = float(info_dct['max_aaf_all'][0])

            if max_aaf_all > max_maf:
                # Variant is frequent
                benign_clinvar = False

                if 'clinvar_sig' in info_dct.keys():
                    # Variant is frequent and in clinvar
                    path_label = "".join(info_dct['clinvar_sig']).lower()
                    if any(path in path_label for path in non_path_lst):
                        # Variant is frequent and not pathogenic in clinvar
                        benign_clinvar = True
                    else:
                        # Variant is frequent but pathogenic in clinvar
                        out_vcf.write(rec.__str__())
                        continue

                if 'cosmic_ids' in info_dct.keys() and not benign_clinvar:
                    # Variant is frequent but found in COSMIC database
                    out_vcf.write(rec.__str__())
                    continue

                else:
                    # Variant is frequent and not in COSMIC or pathogenic in clinvar
                    continue

            else:
                # Variant is infrequent
                out_vcf.write(rec.__str__())

        else:
            # Variant not found in databases
            out_vcf.write(rec.__str__())

    out_vcf.close()

def var_in_loh(loh_df,chrom,start,stop):

    in_loh = False
    for i in loh_df.index:
        reg = loh_df.loc[i]
        if reg['chromosome'] == chrom:
            if start >= reg['start'] and stop <= reg['end']:
                break
    return in_loh
