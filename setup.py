from setuptools import setup, find_packages

setup(
    name="cube toolkit",
    version="1.0",
    description='A tool to calculate CUB indices and optimize sequences',
    url='https://github.com/dyyvgug/CUBE',
    author='Yingying Dong',
    author_email='dyyvgug@gmail.com',
    packages=find_packages(),
    install_requires=['scipy','rpy2','codonw'],
    data_files=[('resource', ['resource/RSCU/*.txt', 'resource/weight/*'])],
)
