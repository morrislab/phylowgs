from setuptools import setup, find_packages

setup(
    name='phylowgs',
    packages=find_packages(),
    description='Application for inferring subclonal composition and evolution from whole-genome sequencing data.',
    keywords=[],
    classifiers=[],
    entry_points={
        'console_scripts': [
            'phylowgs_evolve = phylowgs.evolve:main',
            'phylowgs_write_results = phylowgs.write_results:main',
            'phylowgs_create_inputs = phylowgs.parser.create_phylowgs_inputs:main',
            'phylowgs_parse_cnvs = phylowgs.parser.parse_cnvs:main',
        ],
    },
)

