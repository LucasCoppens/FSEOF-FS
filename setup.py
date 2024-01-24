from setuptools import setup, find_packages

# Read requirements.txt
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='fvseof',
    version='0.1',
    packages=find_packages(),
    install_requires=required,
)