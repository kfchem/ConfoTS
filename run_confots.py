#!/usr/bin/env python3

import subprocess
import sys

steps = [
    'confots_step0.py',
    'confots_step1.py',
    'confots_step2.py',
    'confots_step3.py',
    'confots_step4.py',
    'confots_step5.py',
]

for step in steps:
    print(f'\n=== Running {step} ===')
    result = subprocess.run(['python', step])
    if result.returncode != 0:
        print(f'❌ {step} failed. Stopping.')
        sys.exit(result.returncode)

print('\n✅ All steps completed successfully.')