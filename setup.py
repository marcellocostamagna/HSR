from setuptools import setup, find_packages
import os

# Function to read the __version__ variable from HSR/version.py
def read_version():
    version_file = os.path.join(os.path.dirname(__file__), 'HSR', 'version.py')
    namespace = {}
    with open(version_file, 'r', encoding='utf-8') as f:
        exec(f.read(), namespace)
    return namespace['__version__']

# Use the read_version function to get the version
__version__ = read_version()

with open(os.path.join(os.path.dirname(__file__), 'README.md'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

setup(
    name="HSR",
    version=__version__,
    author="Marcello Costamagna", 
    license="GNU Affero General Public License", 
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
