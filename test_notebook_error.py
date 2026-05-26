import matplotlib.pyplot as plt
from upsetplot import from_contents, plot

site_significant_genes = {'liver': ['A', 'B'], 'axilla': ['B', 'C']}
fig = plt.figure(figsize=(10, 6))
upset_data = from_contents(site_significant_genes)
axes = plot(upset_data, show_counts='%d', sort_by='cardinality', fig=fig)
plt.suptitle("Metastatic Convergence: Upregulated Targets across Niches", fontsize=14, y=1.05)
import io
buf = io.BytesIO()
fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
print("SUCCESS")
