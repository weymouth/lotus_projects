#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
import os
import shutil

def run(n_proc=0,run_folder='test',read_folder=None,test=False):
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
            subprocess.call('trash '+run_folder+'/*', shell=True)
    else:
        print('Creating '+run_folder)
        os.makedirs(run_folder)

    print('Setting up in '+run_folder)
    for file_name in os.listdir('postproc'):
        full_file_name = os.path.join('postproc', file_name)
        if os.path.isfile(full_file_name):
            shutil.copy(full_file_name, run_folder)
    if os.path.isfile('lotus.f90'):
        shutil.copy('lotus.f90', run_folder)

    print('Making executable ')
    os.chdir(run_folder)
    subprocess.call('make -C $MGLHOME/geom_lib/ libgeom.a', shell=True)
    subprocess.call('make -C $MGLHOME/oop/ libfluid.a', shell=True)
    subprocess.call('make -f $MGLHOME/oop/Makefile lotus', shell=True)
    print('Finished executable ')

    if(test):
        print('Finished test ')
        print('Popping back up')
        os.chdir('../.')
        return

    print('Running executable ')
    if read_folder:
        if n_proc==0:
            subprocess.call('time ./lotus '+read_folder, shell=True)
        else:
            subprocess.call('mpirun -n '+str(n_proc)+' ./lotus '+read_folder, shell=True)
    else:
        if n_proc==0:
            subprocess.call('time ./lotus', shell=True)
        else:
            subprocess.call('mpirun -n '+str(n_proc)+' ./lotus', shell=True)

    print('Run all python files for postprocessing')
    for file_name in os.listdir():
        if file_name.endswith('.py'):
            print(file_name)
            subprocess.call('python3 '+file_name, shell=True)

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
