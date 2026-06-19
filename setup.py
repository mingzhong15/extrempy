from setuptools import setup, find_packages
import os

# Read requirements from requirements.txt
def read_requirements():
    """Read requirements from requirements.txt file."""
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    with open(requirements_path, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='extrempy',
    version='0.2.1',
    description='post-processing library for atomistic modeling',
    url='https://github.com/mingzhong15/extrempy/',
    author='Mingzhong',
    author_email='zengqiyu15@163.com',
    license='LGPL',
    packages=find_packages(),
    package_data={
        'extrempy.campaign': ['templates/*.j2'],
    },
    include_package_data=True,
    install_requires=read_requirements(),
    zip_safe=False
)
