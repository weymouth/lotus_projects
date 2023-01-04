import os
import numpy as np

def copy(k_lam):
    os.system(f'mkdir -p {k_lam:.0f}')
    os.system(f'mkdir -p {k_lam:.0f}/smooth')
    os.system(f'mkdir -p {k_lam:.0f}/rough')
    os.system(f'cp changef90.py subOpt lotus.py run.py lotus.f90 {k_lam:.0f}/')


def change_run_py(k_lam):
    with open(f"{k_lam:.0f}/run.py","r") as fileSource:
        fileLines = fileSource.readlines()
    
    fileLines[46] = f"    k_lam = {k_lam}\n"
    
    with open(f"{k_lam:.0f}/run.py","w") as fileOutput:
        fileOutput.writelines(fileLines)


def change_run2d_py(k_lam):
    with open(f"{k_lam:.0f}/run.py","r") as fileSource:
        fileLines = fileSource.readlines()
    
    fileLines[21] = f"    k_lam = {k_lam}\n"
    
    with open(f"{k_lam:.0f}/run.py","w") as fileOutput:
        fileOutput.writelines(fileLines)


def change_subOpt(k_lam):
    with open(f"{k_lam:.0f}/subOpt","r") as fileSource:
        fileLines = fileSource.readlines()
    
    fileLines[5] = f"#SBATCH --job-name={k_lam:.0f}\n"
    
    with open(f"{k_lam:.0f}/subOpt","w") as fileOutput:
        fileOutput.writelines(fileLines)


if __name__ == "__main__":
    k_lams = np.arange(4, 68, 4)
    for k_lam in k_lams:
        copy(k_lam)
        change_run_py(k_lam)
        change_run2d_py(k_lam)
        change_subOpt(k_lam)

        os.chdir(f'{k_lam:.0f}')
        os.system(f'sbatch subOpt')
        os.chdir('./..')