#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
import os
import os,signal
import signal
import shutil

def run(n_proc=0,run_folder='test',read_folder=None):
    "setup and run Lotus using lotus.f90 and the files in postproc"

    print('Number of proccessors :{}'.format(n_proc))
    print('Run folder            :{}'.format(run_folder))

    if read_folder:
        print('Read folder           :{}'.format(read_folder))
    else:
        print('No read folder')

    if os.path.exists(run_folder):
        print('Folder '+run_folder+' exists!')
        if read_folder=='./':
            print('Resuming in place')
        else:
            print('Moving contents to trash')
            subprocess.call('rm -r '+run_folder+'/*', shell=True)
    else:
        print('Creating '+run_folder)
        os.makedirs(run_folder)

    print('Setting up in '+run_folder)
    try:
        for file_name in os.listdir('postproc'):
            full_file_name = os.path.join('postproc', file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, run_folder)
    except FileNotFoundError:
        print("No postprocessing set up")

    if os.path.isfile('lotus.f90'):
        shutil.copy('lotus.f90', run_folder)
        try:
            shutil.copy('converged.py', run_folder)
        except FileNotFoundError:
            print("No stopping criteria, the simulation will run its course")

    print('Making executable ')
    os.chdir(run_folder)
    subprocess.call('make -C $MGLHOME/src/geom/ libgeom.a', shell=True)
    subprocess.call('make -C $MGLHOME/src/oop/ libfluid.a', shell=True)
    subprocess.call('make -f $MGLHOME/src/oop/Makefile lotus', shell=True)
    # subprocess.call('python3 $MGLHOME/bin/exit.py &', shell=True)


    print('Finished executable ')

    print('Running executable ')
    if read_folder:
        if n_proc==0:
            subprocess.run('time ./lotus '+read_folder, shell=True)
        else:
            subprocess.run('mpirun -n '+str(n_proc)+' ./lotus '+read_folder, shell=True)
    else:
        if n_proc==0:
            subprocess.run('time ./lotus', shell=True)
        else:
            subprocess.run('mpirun -n '+str(n_proc)+' ./lotus', shell=True)

    print('Run all python files for postprocessing')
    try:
        for file_name in os.listdir('.'):
            if file_name.endswith('.py') and file_name!='converged.py':
                print(file_name)
                subprocess.call('python3 '+file_name, shell=True)
    except FileNotFoundError:
        print("No postproc folder")

    print('Popping back up')
    os.chdir('../.')


def replace(template,dic):
    """
    Write lotus.f90 file by replacing the dic on the template
    and return a potential folder name
    """
    f1 = open(template,'r')
    f2 = open('lotus.f90','w')
    for text in f1:
        for i, j in dic.items():
            text = text.replace(i, j)
        f2.write(text)
    f1.close()
    f2.close()
    return '_'.join([i+'_'+j for i,j in dic.items()])

