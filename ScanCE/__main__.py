"""Main module."""
import sys
import os
import subprocess

__version__ = 'v3.0'

COMMANDS = {
    'scan_ce':         ('ScanCE_v3.py',        'Unified CE detection for SR/LR/SC (v3, recommended)'),
    'Scan_ce_loose':   ('Scan_ce_loose.py',   'Short-read CE detection (loose, v1, legacy)'),
    'Scan_ce_strict':  ('Scan_ce_strict.py',  'Short-read CE detection (strict, v1, legacy)'),
    'Scan_ce_lr_loose':('Scan_ce_longread_loose.py', 'Long-read CE detection (loose, v1, legacy)'),
    'Scan_ce_lr_strict':('Scan_ce_longread_strict.py','Long-read CE detection (strict, v1, legacy)'),
}

def main():
    file_abs_path = os.path.abspath(os.path.dirname(__file__))

    task = sys.argv[1] if len(sys.argv) > 1 else ''

    if task in COMMANDS:
        script, _ = COMMANDS[task]
        args = sys.argv[2:]
        subprocess.run(['python', f'{file_abs_path}/{script}'] + args, check=True)
    else:
        if task and task not in ('-v', '--version'):
            print(f'\nERROR: Unknown command: {task}')
        print(f'\nProgram:\tScanCE')
        print(f'Version:\t{__version__}')
        print(f'Usage:\t\tScanCE <command> [options]')
        print(f'\nCommands:')
        for cmd, (_, desc) in COMMANDS.items():
            print(f'\t{cmd:<22}{desc}')
        print()


if __name__ == '__main__':
    main()
