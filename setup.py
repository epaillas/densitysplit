"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

setup(
    name='densitysplit',
    version='0.0.1',
    description='icharacterization of density-dependent galaxy clustering',
    url='https://github.com/epaillas/densitysplit', 
    author='Enrique Paillas',
        author_email='enrique.paillas@uwaterloo.ca',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.6, <4',
    install_requires=[
        'numpy',
        'nbodykit'
    ],
    scripts=['bin/densitysplit']
)
