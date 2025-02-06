from setuptools import setup, find_packages

setup(
    name='chemax',
    version='0.1.1',
    packages=find_packages(),  # This finds all packages inside chemax/
    description='A simple Python package for routine echem data analyses',
    author='Max Meade',
    url='https://github.com/meadem/chemax',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    install_requires=[]  # this has caused problems, putting off implementation for now
)
