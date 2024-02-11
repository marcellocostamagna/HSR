from setuptools import setup, find_packages
from hsr.version import __version__
import os

with open(os.path.join(os.path.dirname(__file__), 'README.md'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

setup(
    name="HSR",
    version=__version__,
    long_description=long_description,
    long_description_content_type='text/markdown', 
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "rdkit"
    ],
    entry_points={
        'console_scripts': [
            'hsr = hsr.hsr_cli:main']
    },
)