"""Replot retained panels from bundled data, without genomic downloads."""
from pathlib import Path
import argparse
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'project/manuscript'))
import retained_panels
import build_delta292_helper_abc_figure as helper
import build_delta292_mechanism_discrimination_figure as mechanisms
import build_delta292_why_figure as lesions


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/retained_panels')
    args=parser.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    data=ROOT/'Supplementary_Data'
    retained_panels.draw_eightq_network(data/'Figure_7_solo_LTR_haplotype_counts.tsv',out/'Figure_7.png')
    retained_panels.draw_sevenp_donors(data/'Figure_S6_donor_categories.tsv',out/'Sevenp_donor_categories.png')
    retained_panels.draw_type_state_counts(data/'Table_S11_direct_TypeI_locus_calls.tsv',out/'Figure_S11.png')
    for module,name in [(helper,'Figure_S15'),(mechanisms,'Figure_S14'),(lesions,'Type_I_lesion_context')]:
        module.OUT=out
        module.PNG=out/(name+'.png')
        module.PDF=out/(name+'.pdf')
        module.main()
    print(out)


if __name__=='__main__':
    main()
