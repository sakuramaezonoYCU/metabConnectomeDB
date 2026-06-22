import sys
from IPython.display import Image

img = Image(filename="foo.png") if False else Image(data=b"dummy")
print("hasattr _repr_png_:", hasattr(img, "_repr_png_"))
print("hasattr data:", hasattr(img, "data"))
print("type of data:", type(img.data))
print("_repr_png_ callable?:", callable(getattr(img, "_repr_png_", None)))
print("returns:", type(img._repr_png_()))
