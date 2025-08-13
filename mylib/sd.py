import re
import os
import subprocess

# perlbrew switch perl-5.6.2

# FJO
ASD='UCUUUCCUCCACA'
FREE2BIND_PATH='/fs/ess/PAS2967/dengyw144/Bacillus_subtilis/free2bind'
def run_fs(seq):
    WORKING_DIR=os.getcwd()
    os.chdir(FREE2BIND_PATH)
    with open(f'tmp.fasta','w') as f:
        f.write(f">\n{seq}")
    result=subprocess.run(["perl",
    'free_scan.pl',
    "-O",
    ASD,
    f'tmp.fasta'
    ], capture_output=True, text=True)
    output=result.stdout
    os.chdir(WORKING_DIR)
    return parse_fs(output)


def parse_fs(data):
    '''
    Parameters:
        Raw output string from free_scan
    Return:
        E_min, spacing
    '''
    result = []

    # Regular expression to capture: first float (with optional minus and decimals) and last int
    pattern = re.compile(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\D*#.*?(\d+)\s*$')

    for line in data.strip().splitlines():
        match = pattern.search(line)
        if match:
            num = float(match.group(1))
            idx = int(match.group(2))
            result.append((num, idx))
    try:
        min_tuple = min(result, key=lambda x: x[0])
    except:
        return 0,0
    E_min, spacing=min_tuple
    return E_min, spacing

if __name__=='__main__':
    data = """
    # Leader Length: 0
    # Some header or comment
    0       # A 0
    0       # A 1
    -1.247805       # U 7
    -8.984685       # A 18
    -0.710585000000002      # G 20
    Random text -3.60793        # G 21
    """
    print(parse_fs(data))