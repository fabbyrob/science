from setuptools import setup
import sys

requires = []

# install ordereddict for python < 2.7
py_version = sys.version_info
if py_version.major == 2 and py_version.minor < 7:
    requires.append('ordereddict')

setup(
    name='PyVCF',
    py_modules=['vcf'],
    scripts=['vcf_melt'],
    author='James Casbon',
    author_email='casbon@gmail.com',
    description='Variant Call Format (VCF) parser for python',
    test_suite='test.test_vcf.suite',
    requires=requires
)
