# -*- coding: utf-8 -*-

from setuptools import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
    cmdclass = {'build_ext': build_ext}
except ImportError:
    cmdclass = {}

import numpy

# TODO: how to make this non-install-breaking?
np_include = numpy.get_include()
version = '0.0.2'

with open('README.rst') as src:
    doc = src.read()

with open('requirements.txt') as src:
    requirements = [line.strip() for line in src]

with open('requirements-dev.txt') as src:
    test_requirements = [line.strip() for line in src]

setup(name='align',
      version=version,
      cmdclass=cmdclass,
      description="polite, proper sequence alignment",
      long_description=doc,
      ext_modules=[
          Extension("align/calign",
                    sources=["align/calign.c"],
                    include_dirs=[np_include, "align"])],
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
      install_requires=requirements,
      test_require=test_requirements,
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
