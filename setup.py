from setuptools import setup

setup(
    name="umi-collapser",
    author="Nikolas Barkas",
    author_email="nbarkas@broadinstitute.org",
    description="scRNA-seq UMI read collapser",
    url="https://github.com/broadinstitute/umi-collapser-internal",
    install_requires=[
        "pysam",
        "tqdm"
    ],
    entry_points={
        'console_scripts': [
            'umi-collapser = umi-collapser.__main__:main'
        ]
    }

)