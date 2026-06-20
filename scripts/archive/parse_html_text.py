import os
import glob
from html.parser import HTMLParser

class TextExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_style_or_script = False
        self.text_content = []

    def handle_starttag(self, tag, attrs):
        if tag in ['style', 'script']:
            self.in_style_or_script = True

    def handle_endtag(self, tag):
        if tag in ['style', 'script']:
            self.in_style_or_script = False

    def handle_data(self, data):
        if not self.in_style_or_script:
            stripped = data.strip()
            if stripped:
                self.text_content.append(stripped)

    def get_text(self):
        return '\n'.join(self.text_content)

workspace_dir = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB'
output_dir = os.path.join(workspace_dir, 'output')
html_files = glob.glob(os.path.join(output_dir, '*.html'))

out_file = os.path.join(output_dir, 'parsed_html_text.txt')

with open(out_file, 'w', encoding='utf-8') as out_f:
    for f in html_files:
        if not os.path.exists(f):
            continue
            
        out_f.write(f"========== TEXT FOR {os.path.basename(f)} ==========\n")
        try:
            with open(f, 'r', encoding='utf-8') as inf:
                content = inf.read()
            extractor = TextExtractor()
            extractor.feed(content)
            
            out_f.write(extractor.get_text())
            out_f.write("\n\n")
        except Exception as e:
            out_f.write(f"Error parsing: {e}\n\n")

print(f"Done parsing HTML to {out_file}.")
