from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='chook', 
    version='0.1.0', 
    description='A comprehensive suite for generating binary optimization problems with planted solutions',
    long_description=readme(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    author='Dilina Perera, Inimfon Akpabio',
    author_email='dilinanp@gmail.com',
    python_requires='>=3.4',
    install_requires=['scipy', 'more-itertools'],
    data_files=[('', ['params.in'])],
    entry_points={
        'console_scripts':[
            'chook=chook.__main__:main',
        ] 
    }
    ) 
