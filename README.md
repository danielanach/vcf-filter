# Small toolkit for VCF filtering and small VCF tasks

Dependencies:
* pandas
* pysam

## Examples:

### Counting variants

```
import count_variants

vcf_list = ['sample_1.vcf','sample_2.vcf','sample_3.vcf.gz']
sample_list = ['1','2','3']

df = count_variants(vcf_list, sample_list)

df
```
| sample | snvs | indels |
| --- | --- | --- |
| 1 | 4500 | 300 |
| 2 | 1245 | 23 |
| 3 | 9823 | 133 |


## Post-bcbio somatic prioritization germline filtering

This one is for filtering out likely germline variants that remain post bcbio tumor-only variant calling. Assumes you used this pipeline with the gemini prioritization. It is not ready for use, just a mock up.

**Basic principles are:**

*Population database filtering:*

Filter out variants that are higher than maximum minor allelic (`max_maf`) fraction in population databases annotated with gemini

UNLESS:
* In cosmic and not benign in clinvar
* In cosmic and not in clinvar
* Not in cosmic but not benign in clinvar

*VAF based filtering:*

Filter out variants that are at very high allelic fraction (suggested `max_af` >0.9) AND not in region of LOH. Requires dataframe containing LOH regions defined. It is quite unlikely that a somatic variant will present in every single cell AND on both alleles, thus it is more likely to be a germline variant.. Unless it is in a region of LOH, in which case it is much more probable. 

```
filter_bcbio_somatic(in_vcf_name,
                         out_vcf_name,
                         max_maf,
                         max_af=0.9,
                         loh_sites=[],
                         non_path_lst=['benign'],
                         loh_df=None)
```
