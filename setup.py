#!/usr/bin/env python
__usage__ = "setup.py command [--options]"
__description__ = "standard install script"
__author__ = "Kunal Mehta (kunal.mehta@mail.utoronto.ca)"

#-------------------------------------------------

from setuptools import (setup, find_packages)
import glob

setup(
    name = 'ns-struc-dm',
    version = '0.0',
    url = 'https://github.com/kunalkishanmehta/ns-struc-dm',
    author = __author__,
    author_email = 'kunal.mehta@mail.utoronto.ca',
    description = __description__,
    license = '',
    scripts = glob.glob('bin/*'),
    packages = find_packages(),
    data_files = [],
    requires = [],
)
