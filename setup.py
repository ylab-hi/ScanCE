#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup(
    author="Xiaotong Lu",
    author_email='lu000016@umn.edu',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A computational workflow for cryptic exon identification in RNA-seq data.",
    install_requires=["pysam", "gffutils", "pandas", "biopython","numpy","scipy.stats"],
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
    url='https://github.com/lu000016/ScanCE',
    version='1.0',
    zip_safe=False,
)
