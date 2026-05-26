import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib
matplotlib.use('Agg')
from upsetplot import from_contents, plot
import io

original_text_init = matplotlib.text.Text.__init__
original_set_position = matplotlib.text.Text.set_position

def patched_text_init(self, x=0, y=0, text='', **kwargs):
    try: x = float(x)
    except Exception:
        try: x = x.item()
        except Exception: pass
    try: y = float(y)
    except Exception:
        try: y = y.item()
        except Exception: pass
    original_text_init(self, x=x, y=y, text=text, **kwargs)

def patched_set_position(self, pos):
    x, y = pos
    try: x = float(x)
    except Exception:
        try: x = x.item()
        except Exception: pass
    try: y = float(y)
    except Exception:
        try: y = y.item()
        except Exception: pass
    original_set_position(self, (x, y))

matplotlib.text.Text.__init__ = patched_text_init
matplotlib.text.Text.set_position = patched_set_position

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
    print("Success with Monkey Patch!")
except Exception as e:
    import traceback
    traceback.print_exc()
