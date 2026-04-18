import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# -------------------------
# STEP 0: Create folders
# -------------------------
os.makedirs("../results/plots", exist_ok=True)
os.makedirs("../results/tables", exist_ok=True)

# -------------------------
# STEP 1: Load Mutation Data
# -------------------------
mut_file = "data_mutations.txt"
df = pd.read_csv(mut_file, sep="\t", comment="#")

print("Mutation Data Loaded")
print(df.head())

# -------------------------
# STEP 2: Filter Mutations
# -------------------------
df_filtered = df[df["Variant_Classification"] != "Silent"]

# -------------------------
# STEP 3: Create Mutation Matrix
# -------------------------
mutation_matrix = pd.crosstab(
    df_filtered["Hugo_Symbol"],
    df_filtered["Tumor_Sample_Barcode"]
)

mutation_matrix = (mutation_matrix > 0).astype(int)

print("Mutation matrix shape:", mutation_matrix.shape)

# -------------------------
# STEP 4: Heatmap (Top Genes)
# -------------------------
top_genes = mutation_matrix.sum(axis=1).sort_values(ascending=False).head(20)
subset = mutation_matrix.loc[top_genes.index]

plt.figure(figsize=(12, 8))
sns.heatmap(subset, cmap="Reds")
plt.title("Top Mutated Genes Heatmap")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.savefig("../results/plots/mutation_heatmap.png")
plt.close()

# -------------------------
# STEP 5: Load RNA Data
# -------------------------
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

# -------------------------
# STEP 8: TP53 Analysis
# -------------------------
gene = "TP53"

if gene in mutation_matrix.index:

    mutated_samples = mutation_matrix.columns[mutation_matrix.loc[gene] == 1]
    common_samples = list(set(mutated_samples) & set(rna.columns))

    print("TP53 mutated samples:", len(mutated_samples))
    print("Common samples for TP53:", len(common_samples))

    if len(common_samples) == 0:
        print("No matching samples found. Check IDs.")
    else:
        mut_group = rna[common_samples]
        non_mut_group = rna.drop(columns=common_samples)

        # Mean expression difference
        mut_mean = mut_group.mean(axis=1)
        non_mut_mean = non_mut_group.mean(axis=1)

        diff = mut_mean - non_mut_mean

        top_diff = diff.sort_values(ascending=False).head(10)

        print("\nTop Differential Genes:\n", top_diff)

        # Save results
        top_diff.to_csv("../results/tables/top_diff_genes.csv")

        # Plot
        plt.figure(figsize=(10, 5))
        top_diff.plot(kind="bar")
        plt.title(f"Top Genes Affected by {gene} Mutation")
        plt.xlabel("Genes")
        plt.ylabel("Expression Difference")
        plt.tight_layout()
        plt.savefig("../results/plots/top_diff_genes.png")
        plt.close()

else:
    print(f"{gene} not found in mutation data")

# -------------------------
# STEP 9: Clustermap
# -------------------------
sns.clustermap(subset, cmap="Reds", figsize=(10, 8))
plt.savefig("../results/plots/mutation_clustermap.png")
plt.close()

print("\n✅ Analysis Completed Successfully!")
