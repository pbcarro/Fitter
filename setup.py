
import os
from shutil import copy
from subprocess import run, PIPE
from warnings import warn
from setuptools import setup, find_packages
from setuptools.command.install import install
from distutils.spawn import find_executable

class InstallClass(install):
    def build_static(self):
        # Find if gcc is available on the system
        if find_executable("gcc") is False:
            raise Exception("gcc compiler not found in PATH; library will not be built!")
        else:
            # Copied the compiler flags from Fitter.c
            print("Building static libraries.")
            lib_cmd = "gcc -Wall -o pyfitter/Fitter.so -shared -fPIC -O3 -funroll-loops pyfitter/Fitter.c -lm -lgsl -lgslcblas"
            process = run(lib_cmd.split(), stdout=PIPE, stderr=PIPE)
            if process.returncode != 0:
                warn(f"gcc compilation returned error code {process.returncode}")

    def run(self):
        self.build_static()
        install.run(self)


setup(
    name="pyfitter",
    version="0.1.0",
    description="Python wrapper for fast rigid rotor asymmetric top simulations",
    author="Brandon Carroll",
    packages=find_packages(),
    author_email="paul.carroll@cfa.harvard.edu",
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
    ],
    scripts=["scripts/DatFileBuilder.py", "scripts/EigenValueSolve.py"],
    cmdclass={"install": InstallClass}
)
