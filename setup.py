from setuptools import setup, find_packages

setup(
    name="cube tool",
    version="1.0",
    description='a tool to calculate mCAI and optimize sequences',
    url='https://github.com/dyyvgug/mCAI',
    author='Yingying Dong',
    author_email='dyyvgug@gmail.com',
    packages=find_packages(),
    install_requires=['rpy2'],
    data_files=[('resource', ['resource/RSCU/*.txt', 'resource/weight/*'])],
)
