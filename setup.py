import setuptools

setuptools.setup(
    name="icaparser",
    version="0.2.9",
    url="https://github.com/Bayer-Group/ica-parser",
    author="Henrik Seidel",
    author_email="heseber+github@mailbox.org",
    description="Parser for JSON files from Illumina Connected Annotations analysis pipeline",
    long_description=open("README.md").read(),
    packages=setuptools.find_packages(),
    install_requires=[
        "natsort>=8.1.0",
        "numpy>=1.22.4",
        "pandas>=1.4.2",
        "setuptools>=60.10.0",
        "tqdm>=4.64.0",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
    include_package_data=True,
    package_data={"": ["data/*"]},
)
