from setuptools import setup, find_packages
from asymmetric_uncertainty import __author__, __contact__, __version__

with open("./README.md","r") as f:
    readme = f.read()

with open("./LICENSE","r") as f:
    license = f.read()

if __name__ == '__main__':
    setup(name="asymmetric_uncertainty",
          version=__version__,
          author=__author__,
          author_email=__contact__,
          packages=find_packages(),
          url="https://github.com/cgobat/asymmetric_uncertainty",
          description="A package that implements a class to represent numeric quantities with asymmetric uncertainties.",
          long_description=readme,
          license=license,
          install_requires=["numpy","matplotlib"])