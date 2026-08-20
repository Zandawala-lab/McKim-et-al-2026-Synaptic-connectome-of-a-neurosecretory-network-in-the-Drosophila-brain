# %%
import pandas as pd
import matplotlib
matplotlib.use('Agg')     # Comment out to use plt.show()
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


# %%
# Project folder and path structure ────────────────────────────────────────────
from pathlib import Path

script_dir = Path(__file__).parent          # .../NSC Connectome/scripts/Figure_8
project_root = script_dir.parent.parent     # .../NSC Connectome

input_dir = project_root / "input"
output_dir = project_root / "output" / "Figure_8"

filename = "Figure8_S1C_opto_validation.csv"
input_file = input_dir / filename

output_dir.mkdir(exist_ok=True)
fig_path = output_dir / "Figure_8_S1C" # saved as .pdf

# %%
# Genotype identifiers (substring match)
EXP_KEYWORD  = 'empty'   # control genotype contains this string
# Genotype NOT containing EXP_KEYWORD is treated as experimental

BEHAVIOR_LABELS = [
   'Abd. curl',
   'Abd. curl +\n genital eversion',
   'Genital eversion',
   'Genital Grooming',
   'Ejaculation',
]

# Column suffixes as they appear in the CSV (after fixing typos)
COL_SUFFIXES = ['curl', 'curl+evert', 'evert', 'genital-groom', 'ejac']

STIM_NUMS = [1, 2, 3, 4, 5]
FIGW = 10
FIGH = 8
txtsize = 18
axislabsize = 28
# ─────────────────────────────────────────────────────────────────────────────

df = pd.read_csv(input_file)
# %%
genotypes = df['genotype'].unique().tolist()
ctrl_geno = [g for g in genotypes if EXP_KEYWORD in g][0]
exp_geno  = [g for g in genotypes if EXP_KEYWORD not in g][0]

def get_col(stim, suffix):
    return f'Stim_{stim}_{suffix}'

def sum_behavior(geno):
    sub = df[df['genotype'] == geno]
    return np.array([
        sub[[get_col(s, suf) for s in STIM_NUMS]].values.sum()
        for suf in COL_SUFFIXES
    ], dtype=float)

exp_totals  = sum_behavior(exp_geno)
ctrl_totals = sum_behavior(ctrl_geno)

# %%
fig, ax = plt.subplots(figsize=(FIGW, FIGH))
x = np.arange(len(BEHAVIOR_LABELS))
bar_width = 0.35

ax.bar(x - bar_width / 2, exp_totals,  width=bar_width, color='black', edgecolor='black', linewidth=0.8)
ax.bar(x + bar_width / 2, ctrl_totals, width=bar_width, color='white', edgecolor='black', linewidth=0.8)

ax.set_xticks(x)
#ax.set_xticklabels(BEHAVIOR_LABELS, fontsize=txtsize, rotation=45)
ax.set_xticklabels(BEHAVIOR_LABELS, fontsize=txtsize, rotation=45, ha='right', rotation_mode='anchor')
ax.set_yticklabels([0,25,50,75,100,125,150,175], fontsize=txtsize)

ax.set_ylabel('Observed behaviors', fontsize=axislabsize)
#ax.set_xlabel('Behavior type', fontsize=axislabsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
#plt.show()
# plt.savefig(f'{fig_path}.svg', format='svg', bbox_inches='tight')
plt.savefig(f'{fig_path}.pdf', format='pdf', bbox_inches='tight')
print(f"Saved {fig_path}.pdf")

# %%
