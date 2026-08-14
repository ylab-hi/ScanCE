#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    author="Xiaotong Lu",
    author_email='lu000016@umn.edu',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description="Unified de novo detection and quantification of cryptic exons from short-read, long-read, and single-cell RNA-seq.",
    install_requires=["pysam>=0.19.0", "gffutils>=0.12"],
    entry_points={
        'console_scripts': [
            'ScanCE=ScanCE.__main__:main'
           
        ],
    },
    license="MIT license",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    #keywords='ScanCE',
    name='ScanCE',
    packages=find_packages(include=['ScanCE', 'ScanCE.*']),
    package_data={'ScanCE': ['config.ini', 'config.ini.example', '*.py']},
    url='https://github.com/ylab-hi/ScanCE',
    version='3.0.0',
    zip_safe=False,
)
