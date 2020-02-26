import hail as hl

def calculate_fst(mt, population, maf_threshold=0, compute_average=True):
    # calculate AF in each population
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((mt.variant_qc.AF[0] > maf_threshold) |
                        (mt.variant_qc.AF[1] > maf_threshold),
                        keep=True)

    fst = (mt
           .group_cols_by(mt[population])
           .aggregate(__n_hom_ref=hl.agg.count_where(mt.GT.is_hom_ref()),
                      __n_hom_alt=hl.agg.count_where(mt.GT.is_hom_var()),
                      __n_het_ref=hl.agg.count_where(mt.GT.is_het_ref()),
                      __n_called=hl.agg.count_where(hl.is_defined(mt.GT))))

    fst = fst.annotate_entries(__AF_alt=(2*fst['__n_hom_alt']+fst['__n_het_ref'])/(2*fst['__n_called']),
                               __AF_ref=(2*fst['__n_hom_ref']+fst['__n_het_ref'])/(2*fst['__n_called']))

    # calculate the expected genotypic counts under Hardy-Weinberg Equilibrium
    fst = fst.annotate_entries(__GC_hom_ref_exp=fst['__n_called']*fst['__AF_alt']*fst['__AF_alt'],
                               __GC_hom_alt_exp=2*fst['__n_called']*fst['__AF_alt']*fst['__AF_ref'],
                               __GC_het_ref_exp=fst['__n_called']*fst['__AF_ref']*fst['__AF_ref'])

    # calculate the local observed & expected heterozygosities
    fst = fst.annotate_entries(__H_obs=fst['__n_het_ref']/fst['__n_called'],
                               __H_exp=1-fst['__AF_alt']**2-fst['__AF_ref']**2)

    # calculate the allele frequency over the total population
    # variant_qc.AF = [0.9868421052631579,0.013157894736842105] __AF_ref, __AF_alt

    fst = fst.annotate_rows(H_I=hl.agg.sum(fst['__n_called']*fst['__H_obs'])/fst.variant_qc.n_called,
                            H_S=hl.agg.sum(fst['__n_called']*fst['__H_exp'])/fst.variant_qc.n_called,
                            H_T=1-fst.variant_qc.AF[0]**2-fst.variant_qc.AF[1]**2)
    fst = fst.annotate_rows(F_IS=(fst['H_S']-fst['H_I'])/fst['H_S'],
                            F_ST=(fst['H_T']-fst['H_S'])/fst['H_T'],
                            F_IT=(fst['H_T']-fst['H_I'])/fst['H_T'])

    fst = fst.rows().select("F_IS", "F_ST", "F_IT")

    if compute_average:
        average_fst = fst.aggregate(hl.agg.mean(fst.F_ST))
        return average_fst
    else:
        return fst


