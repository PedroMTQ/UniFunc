import pathlib
from setuptools import setup,find_packages


package_name='UniFunc'
def get_package_version():
    import os
    from sys import platform
    if platform.startswith('win'):
        SPLITTER = '\\'
    else:
        SPLITTER = '/'


    dir_name=os.path.dirname(os.path.abspath(__file__))
    init_path=f'{dir_name}{SPLITTER}{package_name.lower()}{SPLITTER}__init__.py'
    package_version=None
    with open(init_path) as file:
        for line in file:
            if '__version__' in line:
                package_version=line.replace('__version__','')
                package_version=package_version.strip('\n')
                package_version=package_version.strip()
                package_version=package_version.strip('=')
                package_version=package_version.strip()
                package_version=package_version.strip('"')
                package_version=package_version.strip('"')
    return package_version



# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

long_description='UniFunc is a text mining tool that processes and analysis text similarity between a pair of protein function annotations. It is mainly used as a cross-linking mechanism or redundancy elimination tool when processing annotations without any sort of database identifiers.'

setup(
    name=package_name,
    version=get_package_version(),
    author="Pedro Queirós",
    author_email="pdqueiros@gmail.com",
    description="Tool for similarity analysis of protein function annotations.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/PedroMTQ/UniFunc",
    project_urls={
        "Bug Tracker": "https://github.com/PedroMTQ/UniFunc/issues",
    },
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license="MIT",
    include_package_data=True,
    package_data={
        'myapp': ['data/*.txt'],
    },
    install_requires=['nltk','numpy','requests'],
    entry_points={
        "console_scripts": [
            "unifunc=unifunc.__main__:main",
        ],
    },
)
