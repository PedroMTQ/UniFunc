import setuptools

long_description='UniFunc is a text mining tool that processes and analysis text similarity between a pair of protein function annotations. It is mainly used as a cross-linking mechanism or redundancy elimination tool when processing annotations without any sort of database identifiers.'

setuptools.setup(
    name="UniFunc",
    version="1.2",
    author="Pedro Queirós",
    author_email="pdqueiros@gmail.com",
    description="Tool for similarity analysis of protein function annotations.",
    long_description=long_description,
    url="https://github.com/PedroMTQ/UniFunc",
    project_urls={
        "Bug Tracker": "https://github.com/PedroMTQ/UniFunc/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
    install_requires=['nltk','numpy','python>=3.6','requests']
)
