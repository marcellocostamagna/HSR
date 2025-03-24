from setuptools import setup, find_packages

setup(
    name="hsr",
    use_scm_version=True,
    setup_requires=["setuptools>=45", "setuptools_scm>=6.0.1"],
    author="Marcello Costamagna", 
    license="AGPL-3.0",
    description="Hypershape recognition (HSR): a general framework for moment-based similarity measures",
    long_description=open("ReadME.md", encoding="utf-8").read(),
    long_description_content_type='text/markdown',
    packages=find_packages(exclude=["tests", "tests.*"]),
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
