 from distutils.core import setup, Extension
 import distutils.command.bdist_conda

setup(
    name="unifunc",
    version="1.2",
    distclass=distutils.command.bdist_conda.CondaDistribution,
    conda_buildnum=0,
)
