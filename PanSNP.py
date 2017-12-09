# Wrapper function to run pansnp
from .Pansnp import Pansnp_main
import sys

if __name__ == '__main__':
    Pansnp_main.pansnp_main(sys.argv)