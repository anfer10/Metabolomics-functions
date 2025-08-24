# This script is used to build STOCSY traces using raw 1H data 

from pathlib import Path
from STOCSY_1D import stocsy
import os 

# Path where the file is located
folder_path= os.path.dirname(os.path.abspath(__file__)) # Path
folder_path = Path(folder_path)

# STOCSY 1D
stocsy(folder_path)