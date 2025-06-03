from setuptools import setup, find_packages

setup(
    name='MPRAlib',
    version='0.6.5',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    package_data={
        'MPRAlib': ['data/*.txt'],
    },
    install_requires=[
        'numpy',
        'pandas',
        'pysam',
        'click',
        'scikit-learn',
        'scipy',
        'anndata',
        'seaborn',
        'matplotlib',
        'jsonschema',
        'tqdm',
    ],
    extras_require={
        "test": [
            "pytest",
            "coverage",
        ],
    },
    author='Max Schubach',
    author_email='max.schubach@bih-charite.de',
    description='Library to analyze count data of MPRA experiments.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/kircherlab/MPRAlib',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    entry_points={
        "console_scripts": ["mpralib=mpralib.cli:main"]
    },
)
