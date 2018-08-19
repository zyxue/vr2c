#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'pysam>=0.14.1',
    'tqdm>=4.23.4',
    'biopython>=1.72',
    'pandas>=0.23.1',
    'numpy>=1.14.5',
    'matplotlib>=2.2.2',
    'scipy>=1.1.0',
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Zhuyi Xue",
    author_email='zxue@bcgsc.ca',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Visualize read-to-contig alignment",
    entry_points={
        'console_scripts': [
            'vr2c=vr2c.vr2c:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='vr2c',
    name='vr2c',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/zyxue/vr2c',
    version='3.0.0',
    zip_safe=False,
)
