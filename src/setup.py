#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

from setuptools import setup
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()

setup(name='opf_python',
      version='0.0.6',
      description='Optimal kpoint grid generation.',
      long_description=read_md('README.md'),
      author='Wiley S Morgan',
      author_email='wiley.s.morgan@gmail.com',
      url='',
      license='MIT',
      install_requires=[
          "numpy",
          "termcolor",
          "pyparsing",
      ],
      packages=['opf_python'],
      scripts=[],
      classifiers=[
          'Development Status :: 1 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'License :: OSI Approved :: MIT License',          
          'Operating System :: MacOS',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Physics',
      ],
     )
