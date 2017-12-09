# Wrapper function to run pansnp
from Pansnp import Pansnp_main as pm
import sys

if __name__ == '__main__':
    print(sys.path)
    pm.pansnp_main(sys.argv)