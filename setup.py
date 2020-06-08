from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="umi_collapser",
    version="0.0.1",
    author="Nikolas Barkas",
    author_email="nbarkas@broadinstitute.org",
    packages=["umi_collapser"],
    description="scRNA-seq UMI read collapser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/broadinstitute/umi_collapser-internal",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pysam",
        "tqdm",
        "numpy",
        "scipy"
    ],
    entry_points={
        'console_scripts': [
            'umi_collapser = umi_collapser.__main__:main'
        ]
    }
)
