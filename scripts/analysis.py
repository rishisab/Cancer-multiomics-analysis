import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Create output folders
os.makedirs("../results/plots", exist_ok=True)
os.makedirs("../results/tables", exist_ok=True)

# -------------------------
# STEP 1: Load Mutation Data
# -------------------------
mut_file = "../data/mutations.txt"
df = pd.read_csv(mut_file, sep="\t", comment="#")

print("Data Loaded Successfully")
print(df.head())

# -------------------------
# STEP 2: Basic Analysis
# -------------------------
print("\nMutation Types:\n", df["Variant_Classification"].value_counts())
print("\nTop Genes:\n", df["Hugo_Symbol"].value_counts().head(10))

# -------------------------
# STEP 3: Filter Mutations
# -------------------------
df_filtered = df[df["Variant_Classification"] != "Silent"]

# -------------------------
# STEP 4: Mutation Matrix
# -------------------------
mutation_matrix = pd.crosstab(
    df_filtered["Hugo_Symbol"],
    df_filtered["Tumor_Sample_Barcode"]
)

mutation_matrix = (mutation_matrix > 0).astype(int)

# Save matrix
mutation_matrix.to_csv("../results/tables/mutation_matrix.csv")

# -------------------------
# STEP 5: Heatmap
# -------------------------
top_genes = mutation_matrix.sum(axis=1).sort_values(ascending=False).head(20)
subset = mutation_matrix.loc[top_genes.index]

plt.figure(figsize=(10, 8))
sns.heatmap(subset, cmap="Reds")
plt.title("Top Mutated Genes Heatmap")
plt.savefig("../results/plots/mutation_heatmap.png")
plt.close()

# -------------------------
# STEP 6: Load RNA Data
# -------------------------
rna_file = "../data/rna_expression.txt"
rna = pd.read_csv(rna_file, sep="\t", index_col=0)

# Fix sample IDs
rna.columns = [c[:12] for c in rna.columns]

# -------------------------
# STEP 7: Compare TP53
# -------------------------
gene = "TP53"

if gene in mutation_matrix.index:
    mutated_samples = mutation_matrix.columns[mutation_matrix.loc[gene] == 1]

    common_samples = list(set(mutated_samples) & set(rna.columns))

    mut_group = rna[common_samples]
    non_mut_group = rna.drop(columns=common_samples)

    mut_mean = mut_group.mean(axis=1)
    non_mut_mean = non_mut_group.mean(axis=1)

    diff = mut_mean - non_mut_mean

    top_diff = diff.sort_values(ascending=False).head(10)

    # Save results
    top_diff.to_csv("../results/tables/top_diff_genes.csv")

    # Plot
    top_diff.plot(kind="bar", figsize=(10, 5))
    plt.title(f"Top Genes Affected by {gene} Mutation")
    plt.savefig("../results/plots/top_diff_genes.png")
    plt.close()

print("Analysis Completed!")
