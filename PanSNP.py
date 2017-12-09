# Wrapper function to run pansnp
from Pansnp import Pansnp_main as pm
import sys
import os
if __name__ == '__main__':
    print(sys.path)
    sys.path.append(os.path.join("Pansnp"))
    pm.pansnp_main(sys.argv)