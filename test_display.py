from IPython.display import HTML
import json

obj = HTML('<h1>Test</h1>')
print("hasattr _repr_html_:", hasattr(obj, '_repr_html_'))
if hasattr(obj, '_repr_html_'):
    print("_repr_html_():", obj._repr_html_())

print("hasattr data:", hasattr(obj, 'data'))
if hasattr(obj, 'data'):
    print("data:", obj.data)
