#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--srcdir')
parser.add_argument('-i', '--instdir', default='rel')
parser.add_argument('-bt', '--buildtype', default='RelWithDebInfo')
parser.add_argument('-exnet', '--extern-netgen-source', dest='exnet', action='store_true', default=False)

args, rest = parser.parse_known_args()

builddir = os.getcwd()

srcdir = '/home/alex/code/ngsolve/src'
if args.srcdir is not None:
    srcdir += '-' + args.srcdir
srcdir += '/ngsolve'

instdir = '/home/alex/code/ngsolve/inst'
if args.instdir is not None:
    instdir += '-' + args.instdir

outf = open(instdir + '/buildparams.txt', 'w')

cmd = ['cmake']
cmd.append(srcdir)

if args.exnet:
    cmd.append('-DNETGEN_SOURCE_DIR:PATH=' + srcdir + '/netgen')
    cmd.append('-DNETGEN_DIR:PATH=' + instdir)

cmd.append('-DCMAKE_BUILD_TYPE=' + args.buildtype)
cmd.append('-DINSTALL_DIR:PATH=' + instdir)
cmd.append('-DUSE_MKL:BOOL=1')
cmd.append('-DUSE_MUMPS:BOOL=0')
cmd.append('-DUSE_MPI:BOOL=0')
cmd.append('-DUSE_CCACHE:BOOL=1')
cmd.append('-DUSE_UMFPACK:BOOL=1')
cmd.append('-DUMFPACK_INCLUDE_DIR:PATH=/usr/include/suitesparse')
cmd.append('-DUMFPACK_CHOLMOD_LIB:FILEPATH=/usr/lib/x86_64-linux-gnu/libcholmod.so')
cmd.append('-DLIB_MPIF90:FILEPATH=/usr/lib/libmpi_mpifh.so')
cmd.append('-DUMFPACK_CONFIG_LIB:FILEPATH=/usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so')
cmd.append('-DUMFPACK_LIBRARIES:STRING=/usr/lib/x86_64-linux-gnu/libumfpack.so;/usr/lib/x86_64-linux-gnu/libamd.so;/usr/lib/x86_64-linux-gnu/libcolamd.so;/usr/lib/x86_64-linux-gnu/libcholmod.so;/usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so')
cmd.append('-DUMFPACK_LIB:FILEPATH=/usr/lib/x86_64-linux-gnu/libumfpack.so')
cmd.append('-DUMFPACK_COLAMD_LIB:FILEPATH=/usr/lib/x86_64-linux-gnu/libcolamd.so')
cmd.append('-DMKL_ROOT:PATH=/home/alex/intel/mkl')
cmd.append('-DUMFPACK_AMD_LIB:FILEPATH=/usr/lib/x86_64-linux-gnu/libamd.s')
cmd.append('-DVTUNE_INCLUDE_DIR:PATH=/home/alex/intel/vtune_amplifier_xe/include')

outf.write(str(sys.argv)+'\n\n')
outf.write(str(cmd)+'\n\n')
os.chdir(srcdir)
outf.write(str(subprocess.check_output(['git', 'rev-parse', 'HEAD'])))
os.chdir(builddir)
outf.close()

# netgen
# cmd.append('-DUSE_OCC:BOOL=1')


cmd += rest

subprocess.run(cmd)
