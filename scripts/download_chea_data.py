import os
import requests

def download_chea():
    url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ChEA_2022"
    output_dir = os.path.join(os.path.dirname(__file__), "..", "input", "databases", "chea")
    os.makedirs(output_dir, exist_ok=True)
    
    output_path = os.path.join(output_dir, "ChEA_2022.txt")
    
    if os.path.exists(output_path):
        print(f"ChEA 2022 already exists at {output_path}. Skipping download.")
        return output_path
        
    print(f"Downloading ChEA 2022 from {url}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        print(f"Successfully downloaded to {output_path}")
        return output_path
    except Exception as e:
        print(f"Error downloading ChEA 2022: {e}")
        return None
        raise

if __name__ == "__main__":
    download_chea()
