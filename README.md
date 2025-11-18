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

## Running Example Scripts
```bash
python examples/example1.py
python examples/example2.py
```

