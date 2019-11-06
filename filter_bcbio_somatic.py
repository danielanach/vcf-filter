from pysam import VariantFile

def filter_bcbio_somatic(in_vcf_name,
                         out_vcf_name,
                         max_maf,
                         max_af=0.9,
                         loh_df=None):

    vcf_in = VariantFile(in_vcf_name)

    out_vcf = open(out_vcf_name,mode='w')
    out_vcf.write(vcf_in.header.__str__())

    for rec in vcf_in:

        if 'max_aaf_all' in rec.info.keys():
            # Found in population databases

            if is_frequent(rec,max_maf):

                if is_in_clinvar(rec) and not is_benign_in_clinvar(rec):
                    out_vcf.write(rec.__str__())
                    continue

                if is_cosmic_variant(rec,3):
                    out_vcf.write(rec.__str__())
                    continue

                else:
                    continue

        if is_clonal_and_not_loh(rec,max_af,loh_df):
                continue
        else:
            out_vcf.write(rec.__str__())

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
    if 'dbNSFP_clinvar_clnsig' in rec.info.keys():
        return True
    else:
        return False

def is_benign_in_clinvar(rec):

    benign = False

    if '2' in rec.info['dbNSFP_clinvar_clnsig'] or '3' in rec.info['dbNSFP_clinvar_clnsig']:
        benign = True

    return benign

def is_cosmic_variant(rec,min_cnt):

    cosmic = False

    if 'dbNSFP_COSMIC_CNT' in rec.info.keys():
        for c in rec.info['dbNSFP_COSMIC_CNT']:
            if c >= min_cnt:
                cosmic = True
    return cosmic
