import urllib.request
import os

# List of proteins: (ID, Name)
proteins = {
    'P00766': 'Chymotrypsinogen A (Bovine)',
    'P00760': 'Trypsin (Bovine)',
    'P00775': 'Subtilisin Carlsberg (B. licheniformis)',
    'P00974': 'Pancreatic elastase (Porcine)',
    'P15636': 'Thrombin (Human)'
}

if not os.path.exists('data'):
    os.makedirs('data')

all_fasta = ''

for acc, name in proteins.items():
    # URL for obtaining FASTA from UniProt
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.fasta'

    try:
        with urllib.request.urlopen(url) as resp:
            data = resp.read().decode('utf-8')
            all_fasta += data + '\n'
            print(f'Downloaded: {name} ({acc})')
    except Exception as e:
        print(f'Failed to download {acc}: {e}')

# Saving Results
file_path = 'data/serine_proteases.fasta'
with open(file_path, 'w') as f:
    f.write(all_fasta)

print(f'Done! Saved all sequences to {file_path}')