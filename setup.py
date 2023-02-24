from setuptools import setup, find_packages

setup(
    name="sidescan",
    version='0.0.2',
    packages=find_packages(),
    author="Frederico Schmitt Kremer",
    author_email="fred.s.kremer@gmail.com",
    description="Prediction of side effects of drugs",
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    keywords="chemoinformatics",
    entry_points = {'console_scripts':[
        'sidescan-download=sidescan.download:main',
        'sidescan-preprocess=sidescan.preprocess:main',
        'sidescan-train=sidescan.train:main',
        'sidescan-search=sidescan.search:main',
        'sidescan-server=sidescan.server:main'
        ]},
    install_requires = [
        requirement.strip('\n') for requirement in open("requirements.txt")
    ]
)