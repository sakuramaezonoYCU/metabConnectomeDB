import nbformat

def remove_execute():
    nb_path = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/deepdive_23_metabGeneSig.ipynb"
    
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
        
    for cell in nb.cells:
        if cell.cell_type == 'code':
            if 'nbconvert' in cell.source and '--execute' in cell.source:
                cell.source = cell.source.replace('"--execute", ', '')
                cell.source = cell.source.replace("'--execute', ", "")

    with open(nb_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
        
    print("Removed --execute from HTML export cell!")

if __name__ == '__main__':
    remove_execute()
