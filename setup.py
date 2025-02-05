from setuptools import setup, find_packages

setup(
    name='chemax',
    version='0.1.0',
    packages=find_packages(),  # This finds all packages inside chemax/
    install_requires=[],  # List dependencies here if any (what is the proper formatting here?)
    description='A simple Python package for routine echem data analyses',
    author='Max Meade',
    url='https://github.com/meadem/chemax',
    classifiers=[  # These are optional but help others find your package
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Choose license type here
        'Operating System :: OS Independent',
    ],
)
