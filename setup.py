# -*- coding: utf-8 -*-

from setuptools import setup
from distutils.extension import Extension

version = '0.0.2'

with open('README.rst') as src:
    doc = src.read()

with open('HISTORY.rst') as src:
    history = src.read().replace('.. :changelog:', '').strip()

with open('requirements.txt') as src:
    requirements = [line.strip() for line in src]

with open('requirements-dev.txt') as src:
    test_requirements = [line.strip() for line in src]

# Get numpy dependency from requirements
numpy_req, = [req for req in requirements if
              req.lower().startswith("numpy")]

# Attempt to import ~ this should raise an error when doing setup_requires
# (preinstall) but should be ok during actual install
try:
    import numpy
    np_includes = [numpy.get_include()]
except ImportError:
    np_includes = []

# Similar logic to numpy import for Cython

cython_req, = [req for req in requirements
               if req.lower().startswith("cython")]


setup(name='align',
      version=version,
      description='polite, proper sequence alignment',
      long_description=doc + '\n\n' + history,
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
      package_data={'align': ['data/*']},
      setup_requires=["setuptools>=18.0", numpy_req, cython_req],
      install_requires=requirements,
      test_require=test_requirements,
      ext_modules=[
          Extension('align.calign',
                    ['align/calign.pyx'],
                    include_dirs=np_includes)],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering',
          'Topic :: Text Processing'
      ])
