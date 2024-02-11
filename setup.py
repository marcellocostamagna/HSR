from setuptools import setup, find_packages
from hsr.version import __version__

setup(
    name="HSR",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "rdkit"
    ],
    entry_points={
        'console_scripts': [
            'hsr = HSR.hsr_cli:main']
    },

)