# bastress_et_al_2025
Code examples and installation instructions for the method described in Bastress et al. (2025).

## Setup Instructions
### 1. Python environment setup
To keep dependencies isolated, create a virtual environment:
```bash
# Create a virtual environment in the .venv folder
python -m venv .venv

# Activate the virtual environment
# macOS / Linux
source .venv/bin/activate
# Windows (PowerShell)
.venv\Scripts\Activate.ps1

# To deactivate
deactivate
```

### 2. Installing dependencies
Install dependencies from inside the virtual env
```bash
pip install -r requirements.txt
```

### 3. Installing MARS (for example 3)

To compile MARS, follow the [repository instructions](https://github.com/lorrainea/MARS). For convenience, we provide compiled binaries for macOS and Linux in this [forked release](https://github.com/manulera/MARS/releases/tag/v0.2).

## Running Example Scripts

Example 1 and 2 contain basic examples on how to combine pairwise alignments to create a multiple sequence alignment.
```bash
cd examples
python example1.py
python example2.py
python example3.py # outputs example3-MSA.maf
```

