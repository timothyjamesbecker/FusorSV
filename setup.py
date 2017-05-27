from distutils.core import setup
from Cython.Distutils import build_ext
from distutils.extension import Extension
from distutils.version import StrictVersion
from Cython.Build import cythonize



def two_dot(version):
    v = version.split('.')
    return '.'.join(v[0:min(3,len(v))])

# require pysam is pre-installed
try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')
required_pysam_version = '0.9.0'
if StrictVersion(two_dot(pysam.__version__)) < StrictVersion(required_pysam_version) or\
   StrictVersion(two_dot(pysam.__version__)) > StrictVersion(required_pysam_version):
   raise Exception('pysam version == %s is required; found %s' %
                    (required_pysam_version, pysam.__version__))

#require numpy is pre-installed
try:
    import numpy
except ImportError:
    raise Exception('numpy not found; please install numpy first')
required_numpy_version = '1.10.0'
if StrictVersion(two_dot(numpy.__version__)) < StrictVersion(required_numpy_version):
    raise Exception('numpy version >= %s is required; found %s' %
                    (required_numpy_version, numpy.__version__))

# require scipypy is pre-installed
#try:
#    import scipy
#except ImportError:
#    raise Exception('scipy not found; please install scipy first')
#required_scipy_version = '1.16.0'
#if StrictVersion(two_dot(scipy.__version__)) < StrictVersion(required_scipy_version):
#    raise Exception('numpy version >= %s is required; found %s' %
#                    (required_scipy_version, scipy.__version__))    

# require HTSeq is pre-installed
#try:
#    import HTSeq
#except ImportError:
#    raise Exception('HTSeq not found; please install HTSeq first')
#required_HTSeq_version = '0.6.0'
#if StrictVersion(two_dot(HTSeq.__version__)) < StrictVersion(required_HTSeq_version):
#    raise Exception('HTSeq version >= %s is required; found %s' %
#                    (required_HTSeq_version, HTSeq.__version__))  
    
# require bx-python for liftover is pre-installed
#try:
#    import bx
#except ImportError:
#    raise Exception('bx-python not found; please install bx-python first')
#required_bxpython_version = '0.5.0'
#if StrictVersion(two_dot(bx.__version__)) < StrictVersion(required_bxpython_version):
#    raise Exception('bx-python version >= %s is required; found %s' %
#                    (required_bxpython_version, bx.__version__))

# require mygene is pre-installed
#try:
#    import mygene
#except ImportError:
#    raise Exception('mygene not found; please install mygene first')
#required_mygene_version = '3.0.0'
#if StrictVersion(two_dot(mygene.__version__)) < StrictVersion(required_mygene_version):
#    raise Exception('HTSeq version >= %s is required; found %s' %
#                    (required_mygene_version, mygene.__version__)) 

def get_version():
    """Extract version number from source file."""
    from ast import literal_eval
    with open('fusion_utils.pyx') as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("__version__ not found")    
    
cythonize('fusion_utils.pyx')

extensions = [Extension('fusion_utils',
              sources=['fusion_utils.c'],
              libraries=['m'],
              include_dirs=pysam.get_include()+[numpy.get_include()],
              define_macros=pysam.get_defines(),
              extra_compile_args=['-ffast-math'])]

setup(
    name = 'FusorSV',
    version=get_version(),
    author='Timothy Becker',
    author_email='timothyjamesbecker@gmail.com',
    url='https://github.com/timothyjamesbecker/FusorSV',
    license='GPL License',
    description='SV calling data fusion framework',
    classifiers=['Intended Audience :: Developers',
                 'License :: GPL License',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Cython',
                 'Programming Language :: C',
                 'Operating System :: POSIX',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    cmdclass = { 'build_ext': build_ext },
    ext_modules = extensions,
    py_modules = ['FusorSV','fusion_utils','fusor_utils','svu_utils'],
    install_requires = ['pysam>=0.9.0,<0.9.2',
                        'numpy>=0.10.0,<0.11.0',
                        'scipy>=1.16.0,<1.18.0',
                        'HTSeq>=0.6.0<0.7.0',
                        'bx-python>=0.5.0,<0.7.3',
                        'mygene>=3.0.0']
)    
    
