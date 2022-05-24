from setuptools import find_packages, setup
from pathlib import Path

this_directory = Path(__file__).parent
readme = (this_directory / "README.md").read_text()
history = (this_directory / "HISTORY.md").read_text()

setup(
    name="velosearaptor",
    version="0.2.0",
    author="velosearaptor Developers",
    author_email="",
    url="https://github.com/gunnarvoet/velosearaptor/",
    license="",
    # Description
    description="Python functions for interfacing with RDI ADCP data, mostly based on pycurrents",
    long_description=f"{readme}\n\n{history}",
    long_description_content_type='text/markdown',
    # Requirements
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "xarray",
        "matplotlib",
    ],
    extras_require={
    "test": [  # install these with: pip install velosearaptor[test]
        "pytest>=3.8",
    ],
    },
    # Packaging
    packages=find_packages(include=["velosearaptor", "velosearaptor.*"], exclude=["*.tests"]),
    include_package_data=True,
    zip_safe=False,
    platforms=["any"],  # or more specific, e.g. "win32", "cygwin", "osx"
    # Metadata
    project_urls={"Documentation": ""},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
