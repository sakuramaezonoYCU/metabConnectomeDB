import sys
import os
from execute_pancancer_notebooks import execute_and_export

script_dir = os.path.dirname(os.path.abspath(__file__))
visium_nb = os.path.join(script_dir, "visium_spatial_validation.ipynb")
visium_html = os.path.join(os.path.dirname(script_dir), "output", "visium_spatial_validation_report.html")

if os.path.exists(visium_nb):
    execute_and_export(visium_nb, visium_html, "Ovarian Visium Spatial Transcriptomics Validation")
