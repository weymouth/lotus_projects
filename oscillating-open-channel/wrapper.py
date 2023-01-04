import os
import numpy as np

def copy(c):
    os.system(f'mkdir -p {c:.0f}')
    os.system(f'cp changef90.py subOpt lotus.py run.py lotus.f90 save_np_binary.py {c:.0f}/')
    # os.system(f'cp save_np_binary.py {re:.0f}/')


def change_run_py(c):
    with open(f"{c:.0f}/run.py","r") as fileSource:
        fileLines = fileSource.readlines()
    
    fileLines[38] = f"    c = {c}\n"
    
    with open(f"{c:.0f}/run.py","w") as fileOutput:
        fileOutput.writelines(fileLines)


def change_subOpt(c):
    with open(f"{c:.0f}/subOpt","r") as fileSource:
        fileLines = fileSource.readlines()
    
    fileLines[5] = f"#SBATCH --job-name={c:.0f}\n"
    
    with open(f"{c:.0f}/subOpt","w") as fileOutput:
        fileOutput.writelines(fileLines)


if __name__ == "__main__":
    cs = [512, 2048, 4096]
    for c in cs:
        copy(c)
        change_run_py(c)
        # change_subOpt(c)

        # os.chdir(f'{c:.0f}')
        # os.system(f'sbatch subOpt')
        # os.chdir('./..')