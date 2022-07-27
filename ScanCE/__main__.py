"""Main module."""
import sys
import os
import subprocess
import numpy as np
import pandas as pd
__version__ = 'v1.0'
   
def main():
    file_abs_path = os.path.abspath(os.path.dirname(__file__))

    try:
        task = sys.argv[1]
        if task not in ['Scan_ce_loose', 'Scan_ce_strict']:
            if task not in ['-v', '--version']:
                print(f'\nERROR: Unknown usage.')
            raise ValueError
        else:
            if task == 'Scan_ce_loose':
                print('\n')
                args = sys.argv[2:]
                subprocess.run(
                    ['python', f'{file_abs_path}/Scan_ce_loose.py'] + args, check=True)
            elif task == 'Scan_ce_strict':
                print('\n')
                args = sys.argv[2:]
                subprocess.run(
                    ['python', f'{file_abs_path}/Scan_ce_strict.py'] + args, check=True)
    except:
        print(f'\nProgram:\tScanCE')
        print(f'Version:\t{__version__}')
        print(f'Usage:\t\tScanCE <command> [options]')
        print(f'\nCommands:\tScan_ce_loose\t\t\tTool for detecting cryptic exon splicing events in loose standard')
        print(f'\t\tScan_ce_strict\t\tTool for detecting cryptic exon splicing events in strict standard')   
 
            

if __name__ == '__main__':

    main()
    
