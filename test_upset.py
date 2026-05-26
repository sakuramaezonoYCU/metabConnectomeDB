import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from upsetplot import from_contents, plot
import io

site_significant_genes = {
    'liver': ['A', 'B', 'C'],
    'axilla': ['B', 'C', 'D'],
    'chest wall': ['C', 'E']
}

fig = plt.figure(figsize=(10, 6))
upset_data = from_contents(site_significant_genes)
plot(upset_data, show_counts=True, sort_by='cardinality', fig=fig)

buf = io.BytesIO()
try:
    fig.savefig(buf, format='png')
    print("Success")
except Exception as e:
    import traceback
    traceback.print_exc()
