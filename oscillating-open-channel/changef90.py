# -*- coding: utf-8 -*-


def new_f90(S: float, k_lam: int, t_max=10, L=(1, 3, 1), aspect=8):
    with open("lotus.f90","r") as fileSource:
        fileLines = fileSource.readlines()

    fileLines[15] = f"    integer,parameter  :: S = {S}  ! Stokes length\n"
    fileLines[21] = f"    real               :: lambda=S*1./{k_lam}   ! Roughness wavelength\n"

    # Domain test
    fileLines[13] = f"    real               :: m(3) = [{L[0]}*1./{aspect}, 2., {L[2]}*1./{aspect}]\n"
    fileLines[39] = f"    call xg(2)%stretch(n(2),0.,0.,1.*S,{float(L[1])}*S,prnt=root)\n"

    # Cell aspect ratio
    fileLines[38] = f"    xg(1)%h = {float(aspect)}\n"
    fileLines[40] = f"    xg(3)%h = {float(aspect)}\n"

    # Simulation time, number of oscillation cycles
    fileLines[17] = f"    real               :: finish=2*pi*{t_max}*S, print_int=2*pi*0.02*S, init_time=2*pi*3*S\n"
    with open("lotus.f90","w") as fileOutput:
        fileOutput.writelines(fileLines)


