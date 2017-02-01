#!/usr/bin/env python3
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--srcdir")
parser.add_argument("-i", "--instdir", default="rel")
parser.add_argument("-bt", "--buildtype", default="RelWithDebInfo")

args, rest = parser.parse_known_args()

srcdir = "/home/alex/code/ngsolve/src"
if args.srcdir is not None:
    srcdir += "-" + args.srcdir

instdir = "/home/alex/code/ngsolve/inst"
if args.instdir is not None:
    instdir += "-" + args.instdir

cmd = ["cmake"]
cmd.append(srcdir + "/netgen")
cmd.append("-DCMAKE_BUILD_TYPE=" + args.buildtype)
cmd.append("-DINSTALL_DIR:PATH=" + instdir)
cmd.append("-DUSE_JPEG:BOOL=1")
cmd.append("-DUSE_CCACHE:BOOL=1")
cmd.append("-DUSE_MPEG:BOOL=1")
cmd.append("-DUSE_MPI:BOOL=0")
cmd.append("-DUSE_OCC:BOOL=1")

cmd += rest

subprocess.run(cmd)
