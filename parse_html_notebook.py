from bs4 import BeautifulSoup
import sys

with open(sys.argv[1], 'r') as f:
    soup = BeautifulSoup(f, 'html.parser')

# Find the markdown cell containing "3. Machine Learning Models"
headers = soup.find_all(['h2', 'h3', 'h4'])
target_header = None
for h in headers:
    if '3. Machine Learning Models' in h.text:
        target_header = h
        break

if target_header:
    print(f"FOUND HEADER: {target_header.text}")
    current = target_header.find_parent('div', class_='cell')
    if current:
        while current:
            current = current.find_next_sibling('div', class_='cell')
            if current and 'code_cell' in current.get('class', []):
                code = current.find('div', class_='input_area').text
                print("\n--- CODE CELL ---")
                print(code)
                break
    else:
        current = target_header.find_parent('div', class_='jp-Cell')
        while current:
            current = current.find_next_sibling('div', class_='jp-Cell')
            if current and 'jp-CodeCell' in current.get('class', []):
                code = current.find('div', class_='jp-Editor').text
                print("\n--- CODE CELL ---")
                print(code)
                break
