# from distutils.core import setup
# from Cython.Distutils import build_ext
# from distutils.extension import Extension
# import numpy

# setup(
#   ext_modules = [Extension('metro',
#     sources = ['metro.pyx'],
#     libraries=['gsl', 'gslcblas', 'm'],
#     include_dirs=[numpy.get_include()])],
#   cmdclass = {'build_ext' : build_ext}
# )


import numpy as np

import os
import sys
import subprocess

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


cmdclass = {}
ext_modules = []
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True


args = sys.argv[1:]

base_dir = "particle_packing"

# Make a `cleanall` rule to get rid of intermediate and library files
if "cleanall" in args:
    print "Deleting cython files..."
    # Just in case the build directory was created by accident,
    # note that shell=True should be OK here because the command is constant.
    subprocess.Popen("rm -rf build", shell=True, executable="/bin/bash")
    subprocess.Popen("rm -rf " + base_dir + "/cython/*.c", shell=True, executable="/bin/bash")
    subprocess.Popen("rm -rf " + base_dir + "/ext/*.so", shell=True, executable="/bin/bash")

    # Now do a normal clean
    sys.argv[1] = "clean"

# We want to always use build_ext --inplace
if args.count("build_ext") > 0 and args.count("--inplace") == 0:
    sys.argv.insert(sys.argv.index("build_ext")+1, "--inplace")


# Only build for 64-bit target
os.environ['ARCHFLAGS'] = "-arch x86_64"

if use_cython:
    ext = ".pyx"
    cmdclass={'build_ext': build_ext}
else:
    ext = ".c"

sphere_ext = Extension("particle_packing.ext.sphere",
                   [base_dir + "/cython/sphere_ext" + ext],
                   include_dirs=[np.get_include()],
                   libraries=['gsl', 'gslcblas', 'm'],
                   #extra_compile_args=["-g"],
                   #extra_link_args=["-g"]
                   )


circle_ext = Extension("particle_packing.ext.circle",
                   [base_dir + "/cython/circle_ext" + ext],
                   include_dirs=[np.get_include()],
                   libraries=['gsl', 'gslcblas', 'm'],
                   #extra_compile_args=["-g"],
                   #extra_link_args=["-g"]
                   )


ellipse_ext = Extension("particle_packing.ext.ellipse",
                   [base_dir + "/cython/ellipse_ext" + ext],
                   include_dirs=[np.get_include()],
                   libraries=['gsl', 'gslcblas', 'm'],
                   #extra_compile_args=["-g"],
                   #extra_link_args=["-g"]
                   )

boxcar_ext = Extension("particle_packing.ext.boxcar",
                   [base_dir + "/cython/boxcar_ext" + ext],
                   include_dirs=[np.get_include()],
                   libraries=['gsl', 'gslcblas', 'm'],
                   #extra_compile_args=["-g"],
                   #extra_link_args=["-g"]
                   )

ext_modules += [sphere_ext, circle_ext, ellipse_ext, boxcar_ext]



setup(cmdclass=cmdclass,
      ext_modules=ext_modules)
