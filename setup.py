import setuptools, os

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="siman",
    version="1.0.7",
    author="Dmitry Aksenov",
    author_email="dimonaks@gmail.com",
    description="Manager for DFT calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dimonaks/siman",
    license = 'GPL',
    keywords='DFT VASP Density Functional Theory NEB DFT+U PAW GGA Monte-Carlo',
    # package_dir={'':'siman'},
    packages=setuptools.find_packages(),
    # packages=['siman'],
    # data_files=[(os.path.expanduser("~"), ['siman/simanrc.py'])],
    install_requires=['numpy', 'tabulate', 'pymatgen', 'pandas', 'scipy', 'six', 'matplotlib', 'ase', 'paramiko', 'adjustText'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
)
