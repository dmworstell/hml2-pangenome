"""Draw the retained descriptive panels from their checked source tables."""
from pathlib import Path
import csv
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Circle
import numpy as np


def read_tsv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def observed_multicopy_numbers(per_locus):
    return sorted({cn for counts in per_locus.values() for cn in counts if cn >= 2})


def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, facecolor='white')
    fig.savefig(path.with_suffix('.pdf'), facecolor='white')
    plt.close(fig)
    return path


def draw_eightq_network(source_table, path):
    network = read_tsv(source_table)
    if len(network) != 13 or sum(int(r['observations']) for r in network) != 461:
        raise ValueError('Unexpected 8q11.23 sequence denominator or haplotype count')
    if int(network[0]['observations']) != 432:
        raise ValueError('Unexpected dominant 8q11.23 haplotype count')
    fig, (ax, net) = plt.subplots(1, 2, figsize=(10, 4.9),
        gridspec_kw={'width_ratios': [1, 1.15]}, layout='constrained')
    ax.set(xlim=(0, 10), ylim=(0, 10))
    ax.axis('off')
    ax.text(0, 9.8, 'A', fontsize=15, weight='bold')
    ax.text(5, 9.8, '8q11.23 provirus', ha='center', fontsize=15)
    ax.text(5, 8.7, '1 provirus / 583 solo-LTRs', ha='center', fontsize=11)
    ax.plot([.8, 9.2], [7.8, 7.8], color='#999999', lw=6)
    for x in [.8, 8.6]:
        ax.add_patch(Rectangle((x, 7.3), .6, 1, color='#62666A'))
        ax.text(x+.3, 7.8, 'LTR', ha='center', va='center', color='white', fontsize=8)
    for x in [.35, 9.65]:
        ax.text(x, 7.35, 'CACAC', ha='center', va='center', fontsize=8, color='#5E6065')
    for label, x, w, y, color in [
        ('gag', 1.4, 2.7, 6.2, '#507EA1'), ('pro', 3.2, 1.6, 5., '#619A54'),
        ('pol', 3.9, 3.9, 3.8, '#DAB345'), ('Δenv', 6.8, 2., 2.6, '#D56464'),
        ('np9', 6.85, 1.15, 1.4, '#8C729F')]:
        ax.add_patch(Rectangle((x, y), w, .5, color='#62666A'))
        ax.text(x+w/2, y+.25, label, ha='center', va='center', color='white', fontsize=10)
    for x, y, label in [(2.8, 6.7, 'stop'), (5.65, 4.3, 'stop'),
                        (6.15, 3.8, 'FS'), (7.6, 1.4, 'stop')]:
        direction = -1 if label == 'FS' or y < 2 else 1
        ax.plot([x, x], [y, y+direction*.45], color='#934D58', lw=1.5)
        ax.text(x, y+direction*.65, label, ha='center', va='center', fontsize=9, color='#934D58')
    ax.plot([6.7, 6.7], [1.1, 8.2], ls='--', color='#5E6065', lw=1.2)
    ax.text(6.7, 8.35, 'Δ292', ha='center', fontsize=11)
    ax.text(5, .15, 'The two 968-bp LTRs are identical.', ha='center', fontsize=10)
    net.set(xlim=(-1.75, 1.75), ylim=(-1.65, 1.65), aspect='equal')
    net.axis('off')
    net.text(-1.72, 1.61, 'B', fontsize=15, weight='bold')
    net.text(0, 1.61, 'Solo-LTR haplotype network', ha='center', fontsize=15)
    palette = {'AFR': '#0072B2', 'AMR': '#E69F00', 'EAS': '#7A5195',
               'EUR': '#009E73', 'SAS': '#CC79A7', 'Unknown': '#888888'}
    angles = np.linspace(0, 2*np.pi, len(network)-1, endpoint=False)
    for i, row in enumerate(network):
        center = (0, 0) if i == 0 else (1.24*np.cos(angles[i-1]), 1.24*np.sin(angles[i-1]))
        count = int(row['observations'])
        if sum(int(row[pop]) for pop in palette) != count:
            raise ValueError('8q11.23 population counts do not sum to node size')
        radius = .55 if i == 0 else .09 + .06*np.sqrt(count)
        if i:
            net.plot([0, center[0]], [0, center[1]], color='#777777', lw=.8, zorder=0)
            net.text(center[0]*.69, center[1]*.69, str(row['substitutions_from_dominant']),
                ha='center', va='center', fontsize=8,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .2})
        if i == 0:
            # All spokes terminate on the boundary of the whole sequence node.
            net.add_patch(Circle(center, radius + .055, facecolor='white',
                edgecolor='#777777', linewidth=1.1, zorder=1))
        theta = 0
        for pop, color in palette.items():
            theta2 = theta + 360*int(row[pop])/count
            if theta2 > theta:
                net.add_patch(Wedge(center, radius, theta, theta2, facecolor=color,
                    edgecolor='white', linewidth=.3, zorder=2))
            theta = theta2
        net.text(*center, str(count), ha='center', va='center', fontsize=10 if i == 0 else 8,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .5})
    net.legend(handles=[Rectangle((0,0), 1, 1, facecolor=c, label=k) for k,c in palette.items()
        if any(int(row[k]) for row in network)],
        loc='lower center', bbox_to_anchor=(.5, -.07), ncol=3, frameon=False, fontsize=8)
    return save(fig, path)



def draw_type_state_counts(source_table, path):
    rows = sorted((r for r in read_tsv(source_table) if r['locus'] != 'TOTAL'),
                  key=lambda r: int(r['callable_typeI']), reverse=True)
    if len(rows) != 20 or sum(int(r['callable_typeI']) for r in rows) != 9733:
        raise ValueError('Unexpected Type-I denominator')
    if any(int(r['callable_typeII']) != 0 for r in rows):
        raise ValueError('A Type-II call is present, review the figure')
    fig, ax = plt.subplots(figsize=(8.5, 6.1), layout='constrained')
    counts = [int(r['callable_typeI']) for r in rows]
    y = np.arange(len(rows))
    ax.barh(y, counts, color='#163B75')
    ax.set_yticks(y, [r['locus'].removesuffix('_new') for r in rows], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, 630)
    for yy, n in zip(y, counts):
        ax.text(n+5, yy, str(n), va='center', fontsize=9)
    ax.set_xlabel('Callable internal-bearing alleles, all Type I', fontsize=11)
    ax.set_title('Type calls at twenty loci with Type-I alleles', fontsize=15)
    ax.text(.98, .05, 'Retained Type II: 0 at every locus', transform=ax.transAxes, ha='right', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    return save(fig, path)


def draw_sevenp_donors(source_table, path):
    rows = read_tsv(source_table)
    counts = Counter(r['category'] for r in rows)
    if len({r['donor'] for r in rows}) != 292 or sorted(counts.values()) != [7, 26, 259]:
        raise ValueError('Unexpected 7p22.1 donor categories')
    fig, ax = plt.subplots(figsize=(9, 6.1))
    colors = ['#FFC107']*7 + ['#5A9BD4']*259 + ['#6E6E6E']*26
    for i, color in enumerate(colors):
        ax.add_patch(Rectangle((i % 20, i // 20), .96, .96, facecolor=color))
    ax.set(xlim=(-.2, 20), ylim=(-.2, 15), aspect='equal')
    ax.axis('off')
    ax.set_title('7p22.1 coding states in 292 donors', fontsize=16, pad=14)
    labels = ['Four intact ORFs, no Y195C (7/292)',
              'Four intact ORFs, with Y195C (259/292)',
              'No observed copy meeting either criterion (26/292)']
    handles = [Rectangle((0,0), 1, 1, facecolor=c, label=l) for c,l in
               zip(['#FFC107', '#5A9BD4', '#6E6E6E'], labels)]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5,-.01),
              frameon=False, fontsize=11)
    return save(fig, path)


def draw_duplicated_groups(path):
    fig, ax = plt.subplots(figsize=(7, 5.35))
    ax.set(xlim=(0,1), ylim=(0,1))
    ax.axis('off')
    ax.text(.52,.97,'Duplicated HML-2 groups',ha='center',va='top',fontsize=16)
    groups = [('Telomeric Type I', ['13p13','15p13b']),
              ('Telomeric Type II / 4q35.2', ['15p13a','21p13','22p13','4q35.2']),
              ('1p36.21', ['a','b','c']), ('8p23.1', ['b','c','d','e']), ('Xq28', ['a','b'])]
    for i, (title, labels) in enumerate(groups):
        y=.78-i*.165
        ax.text(.08,y+.055,title,fontsize=12)
        x=np.linspace(.15,.9,len(labels))
        ax.plot([x[0],x[-1]],[y,y],color='#737D83',lw=1.2)
        for xx,label in zip(x,labels):
            ax.text(xx,y,label,ha='center',va='center',fontsize=12,
                    bbox={'facecolor':'white','edgecolor':'#737D83','pad':3})
    return save(fig,path)
