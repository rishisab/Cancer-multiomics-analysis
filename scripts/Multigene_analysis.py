# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:54:53 2026

@author: 91907
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

os.makedirs("../results/plots",exist_ok=True)
os.makedirs("../results/tables",exist_ok=True)

mut_file="data_mutations.txt"
df=pd.read_csv(mut_file,sep="\t",comment="#")
print(df.head())

df_filter=df[df["Variant_Classification"]!="Silent"]

mutation_matrix=pd.crosstab(
    df_filter["Hugo_Symbol"],
    df_filter["Tumor_Sample_Barcode"])
mutation_matrix=(mutation_matrix > 0).astype(int)
print("Mutation Matrix Shape",mutation_matrix.shape)
top_genes = mutation_matrix.sum(axis=1).sort_values(ascending=False).head(20)
print(top_genes)
subset = mutation_matrix.loc[top_genes.index]
plt.figure(figsize=(12, 8))
sns.heatmap(subset, cmap="Reds")
plt.title("Top Mutated Genes Heatmap")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.savefig("../results/plots/mutation_heatmap.png")
plt.close()
rna_file = "data_mrna_seq_tpm.txt"
rna = pd.read_csv(rna_file, sep="\t", comment="#", index_col=0)

print("RNA Data Loaded:", rna.shape)

# -------------------------
# STEP 6: Fix Sample IDs
# -------------------------
# Convert both to 12-character TCGA IDs
rna.columns = [c[:12] for c in rna.columns]
mutation_matrix.columns = [c[:12] for c in mutation_matrix.columns]

# Remove duplicates after trimming
rna = rna.groupby(rna.columns, axis=1).mean()
mutation_matrix = mutation_matrix.groupby(mutation_matrix.columns, axis=1).max()

# -------------------------
# STEP 7: Find Common Samples
# -------------------------
common_samples = list(set(rna.columns) & set(mutation_matrix.columns))

print("Mutation samples:", len(mutation_matrix.columns))
print("RNA samples:", len(rna.columns))
print("Common samples:", len(common_samples))

# Save mutation matrix
mutation_matrix.to_csv("../results/tables/mutation_matrix.csv")
all_results=[]
for gene in mutation_matrix.index:

    try:
        mutated_samples = mutation_matrix.columns[mutation_matrix.loc[gene] == 1]

        # Find common samples
        common_samples = list(set(mutated_samples) & set(rna.columns))

        # Skip if no overlap
        if len(common_samples) < 3:
            continue

        mut_group = rna[common_samples]
        non_mut_group = rna.drop(columns=common_samples)

        # Skip if groups are too small
        if mut_group.shape[1] < 2 or non_mut_group.shape[1] < 2:
            continue

        mut_mean = mut_group.mean(axis=1)
        non_mut_mean = non_mut_group.mean(axis=1)

        diff = mut_mean - non_mut_mean

        top_genes = diff.sort_values(ascending=False).head(5)
        print(f"\nGene: {gene}")
        print("Top affected genes:", list(top_genes.index))

        # Store results
        for tg, val in top_genes.items():
            all_results.append({
                "Mutated_Gene": gene,
                "Affected_Gene": tg,
                "Expression_Diff": val
            })

    except Exception as e:
        print(f"Error with {gene}: {e}")
        
results_df = pd.DataFrame(all_results)
results_df.to_csv("../results/tables/all_gene_analysis.csv", index=False)
top_summary = results_df.groupby("Mutated_Gene")["Expression_Diff"].mean().sort_values(ascending=False).head(10)

top_summary.plot(kind="bar", figsize=(10,5))
plt.title("Top Impactful Mutated Genes")
plt.tight_layout()
plt.show()
