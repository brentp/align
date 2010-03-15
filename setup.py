from setuptools import setup, find_packages
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
    cmdclass = {'build_ext': build_ext}
except ImportError:
    cmdclass = None

version = '0.0.1'
import numpy
np_include = numpy.get_include()
doc = open('README.rst').read()

setup(name='align',
      version=version,
      cmdclass=cmdclass,
      description="polite, proper sequence alignment",
      long_description=doc,
      ext_modules=[ Extension("align/calign", sources=["align/calign.c"], include_dirs=[np_include, "align"])],
      keywords='sequence bioinformatics alignment text',
      url='http://github.com/brentp/align/',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='BSD',
      test_suite='nose.collector',
      include_package_data=True,
      zip_safe=False,
      packages=['align'],
      package_dir={'align': 'align'},
      install_requires=['numpy'],
      #entry_points= { 'console_scripts': ['align = align:main'] },
    classifiers   = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering',
        'Topic :: Text Processing'
        ],
)
