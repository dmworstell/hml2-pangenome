"""Rebuild only Figures 1 and 3 and preserve the regional comparison evidence."""
from pathlib import Path
import csv
import runpy
import sys
from Bio import SeqIO
import build_narrative_main_figures as owner

PROJECT=Path(__file__).parent
OUT=PROJECT/'figures/comment_corrections_20260914'
OUT.mkdir(parents=True,exist_ok=True)
owner.OUT=OUT
owner.SUPPLEMENT=OUT/'structural_tables'
owner.SUPPLEMENT.mkdir(exist_ok=True)
rows,roster=owner.load_catalog()
owner.build_figure_1(rows,roster)
runpy.run_path(str(PROJECT/'restore_figure3_20260914.py'),run_name='__main__')

source=PROJECT.parent/'results/acroc_resolved_20260914/phylogeny/selected_sequence_clusters.fasta'
seq={x.id:str(x.seq) for x in SeqIO.parse(source,'fasta')}
result=[]
for region in ['LTR','pol','env']:
    query=next(k for k in seq if k.startswith(region+'|4q35.2_hg38__hap1'))
    for locus in ['15p13a','21p13','22p13']:
        values=[]
        for target in seq:
            if not target.startswith(region+'|'+locus+'__'):continue
            a,b=seq[query],seq[target]
            assert len(a)==len(b)
            n=sum(x in 'ACGT' and y in 'ACGT' for x,y in zip(a,b))
            m=sum(x!=y and x in 'ACGT' and y in 'ACGT' for x,y in zip(a,b))
            gaps=sum((x=='-')!=(y=='-') for x,y in zip(a,b))
            values.append((m/n,m,n,gaps,target))
        distance,m,n,gaps,target=min(values)
        result.append({'region':region,'query':query,'target':target,'mismatches':m,'jointly_called_bases':n,'one_sided_gap_columns':gaps,'identity_at_jointly_called_bases':1-distance})
assert any(x['region']=='pol' and x['target'].startswith('pol|21p13') and x['mismatches']==4 and x['jointly_called_bases']==2741 for x in result)
with (OUT/'regional_sequence_comparisons.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
print('Figures and regional sequence evidence saved')
