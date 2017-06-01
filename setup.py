#!/usr/env/bin/python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import pysam
import numpy

def two_dot(version):
    v = version.split('.')
    return '.'.join(v[0:min(3,len(v))])

def get_version():
    """Extract version number from source file."""
    from ast import literal_eval
    with open('fusorsv/fusion_utils.pyx') as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("__version__ not found")    

cythonize('fusorsv/fusion_utils.pyx')
extensions = [Extension('fusion_utils',
                        sources=['fusorsv/fusion_utils.pyx'],
                        include_dirs=pysam.get_include()+[numpy.get_include()],
                        define_macros=pysam.get_defines(),
                        extra_compile_args=['-ffast-math'])]

setup(
    name = 'fusorsv',
    version=get_version(),
    author='Timothy Becker',
    author_email='timothyjamesbecker@gmail.com',
    url='https://github.com/timothyjamesbecker/FusorSV',
    license='GPL 3 License',
    description='SV calling data fusion framework',
    classifiers=['Intended Audience :: Developers',
                 'License :: GPL 3 License',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Cython',
                 'Programming Language :: C',
                 'Operating System :: POSIX',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    cmdclass = { 'build_ext': build_ext },
    ext_modules = extensions,
    packages =   ['fusorsv'],
    scripts    = ['bin/FusorSV.py'])#,
    #install_requires = [])
                        #'cython>=0.24.0,<0.25.0',
                        #'numpy>=0.10.0,<0.12.0',
                        #'pysam>=0.9.0,<0.9.2',
                        #'bx-python>=0.5.0,<0.7.3', #now optional
                        #'mygene>=3.0.0']           #now optional
