"""Render the retained supplementary trees and mark loci discussed in Results.

No tree inference or sequence analysis is performed. The existing page renderer,
source trees, locus counts and subfamily assignments are retained.
"""
from pathlib import Path
import argparse, csv, hashlib, inspect, json, sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'project/manuscript'))
import build_narrative_main_figures as owner

parser=argparse.ArgumentParser()
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
OUT=args.out.resolve();OUT.mkdir(parents=True,exist_ok=True)
PHY=ROOT/'Supplementary_Data/Phylogeny'
owner.PHY=PHY
owner.TREE_LTR5HS="#007F73"
owner.TREE_LTR5A="#CEAA36"
owner.TREE_LTR5B="#D7733F"
owner.SUBFAMILY_AUTHORITY=Path(__file__).resolve().parent/'inputs/type1/subfamily_authority.tsv'
count_authority=json.loads((PHY/'count_authority.json').read_text())
# Page geometry matches the retained supplemental rendering in
# project/results/acroc_resolved_20260914/phylogeny/rebuild_assignment_phylogenies.py.
renderer=inspect.getsource(owner.build_representative_phylogeny)
renderer=renderer.replace('        "LTR5": PURPLE,','        "LTR5": PURPLE,\n        "non-LTR5": "#343A40",')
renderer=renderer.replace('    x_label = x_tree * 1.07','    x_label = x_tree * 1.03')
renderer=renderer.replace('    for (tip_y, tip_x, label, color), target_y in zip(labels, label_y):','    for label_index, ((tip_y, tip_x, label, color), target_y) in enumerate(zip(labels, label_y)):\n        target_x = x_label + (label_index % 2) * x_tree * 0.33')
renderer=renderer.replace('            [tip_x, x_label * 0.985],\n            [tip_y, target_y],','            [tip_x, x_label * 0.985, target_x * 0.985],\n            [tip_y, target_y, target_y],')
renderer=renderer.replace('            x_label,\n            target_y,','            target_x,\n            target_y,')
renderer=renderer.replace('    ax.set_xlim(-0.065 * x_tree, x_tree * 1.45)','    ax.set_xlim(-0.065 * x_tree, x_tree * 1.73)')
renderer=renderer.replace('    path = PHY / output_name', '''    fig.subplots_adjust(top=0.956, bottom=0.095)
    fig.text(0.51, 0.046, 'Boxed labels mark the loci discussed in the regional tree comparison', ha='center', fontsize=7, color='#665020')
    fig.text(0.51, 0.024, 'n = retained source-bound sequences; up to two modal clusters per locus', ha='center', va='bottom', fontsize=7)
    fig.text(0.51, 0.006, 'Teal: LTR5Hs   Gold: LTR5A   Orange: LTR5B; HML-11 comparators excluded', ha='center', va='bottom', fontsize=7)
    marked=[]
    for text in ax.texts:
        if text.get_text().split(' (n=')[0] in MARKED_LOCI:
            text.set_fontweight('bold')
            text.set_bbox({'boxstyle':'round,pad=0.12','facecolor':'#FFF1CD','edgecolor':'#8A6824','linewidth':0.8})
            marked.append(text.get_text().split(' (n=')[0])
    assert set(marked)==MARKED_LOCI, (output_name, marked, MARKED_LOCI)
    fig.canvas.draw()
    boxes=[(text.get_text(),text.get_window_extent(fig.canvas.get_renderer())) for text in ax.texts if '(n=' in text.get_text()]
    for i,(a,box) in enumerate(boxes):
        for b,other in boxes[i+1:]:
            assert not box.overlaps(other),(output_name,a,b)
    RENDER_QA.append({'panel':output_name,'labels':len(boxes),'marked_loci':sorted(marked),'overlapping_labels':0})
    path = OUTPUT / output_name''')
renderer=renderer.replace('    plt.close(fig)\n    return path','    fig.savefig(path.with_suffix(".svg"), facecolor="white")\n    fig.savefig(path.with_suffix(".pdf"), facecolor="white")\n    plt.close(fig)\n    return path')
owner.OUTPUT=OUT;owner.RENDER_QA=[]
exec(compile(renderer,'retained_supplemental_tree_renderer','exec'),owner.__dict__)
panels=[('gag','Gag','Figure_S7A',{'19p12c','8p23.1a'}),
        ('pro','Pro','Figure_S7B',{'19p12c','7p22.1'}),
        ('LTR','LTR','Figure_S8A',{'19p12c','10q24.2','12q14.1'}),
        ('env','Env','Figure_S8B',{'19p12c','10p12.1'}),
        ('pol','Pol','Figure_S8C',{'19p12c','10q24.2','22q11.21','6q14.1'})]
inputs={}
for region,title,name,marked in panels:
 filename='hml2_pan_ltr_expanded_tree.nwk' if region=='LTR' else f'hml2_pan_orf_{region}_tree.nwk'
 inputs[filename]=hashlib.sha256((PHY/filename).read_bytes()).hexdigest()
 owner.MARKED_LOCI=marked
 owner.build_representative_phylogeny(filename,title,name+'.png',mark_ltr_clades=True,
     count_authority=count_authority[region],figsize=(7.1,9.0),label_fontsize=9.0,
     title_fontsize=11,axis_fontsize=9,clade_label_fontsize=8,star_size=100)
 assert hashlib.sha256((PHY/filename).read_bytes()).hexdigest()==inputs[filename]
(OUT/'supplemental_tree_verification.json').write_text(json.dumps({'status':'PASS','tree_inputs_unchanged':inputs,'rendered_panels':owner.RENDER_QA},indent=2)+'\n')
print(json.dumps(owner.RENDER_QA))
