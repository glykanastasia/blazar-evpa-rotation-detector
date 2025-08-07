"""
Setup script for EVPA Rotation Analysis package.
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    if os.path.exists("README.md"):
        with open("README.md", "r", encoding="utf-8") as fh:
            return fh.read()
    return ""

# Read requirements from requirements.txt
def read_requirements():
    requirements = []
    if os.path.exists("requirements.txt"):
        with open("requirements.txt", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    requirements.append(line)
    return requirements

setup(
    name="evpa-rotation-analysis",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python package for detecting EVPA rotations in polarimetric monitoring data",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/evpa-rotation-analysis",
    
    # Package discovery
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    
    # Dependencies
    install_requires=read_requirements(),
    
    # Optional dependencies
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
        "plotting": [
            "scienceplots>=2.0.0",
        ],
    },
    
    # Entry points for command-line scripts
    entry_points={
        "console_scripts": [
            "evpa-analyze=evpa_rotation.main:main",
        ],
    },
    
    # Package metadata
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    
    python_requires=">=3.8",
    
    # Include additional files
    include_package_data=True,
    package_data={
        "evpa_rotation": ["*.txt", "*.md"],
    },
    
    # Keywords for PyPI search
    keywords="astronomy polarimetry EVPA rotation blazars AGN time-series",
    
    # Project URLs
    project_urls={
        "Bug Reports": "https://github.com/yourusername/evpa-rotation-analysis/issues",
        "Source": "https://github.com/yourusername/evpa-rotation-analysis",
        "Documentation": "https://evpa-rotation-analysis.readthedocs.io/",
    },
)
