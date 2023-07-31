from setuptools import setup, find_packages
import asymmetric_uncertainty.meta as module_metadata

with open("./README.md","r") as f:
    readme = f.read()

with open("./LICENSE","r") as f:
    license = f.read()

if __name__ == '__main__':
    setup(name="asymmetric_uncertainty",
          version=module_metadata.__version__,
          author=module_metadata.__author__,
          author_email=module_metadata.__contact__,
          packages=find_packages(),
          url="https://github.com/cgobat/asymmetric_uncertainty",
          description="A package that implements a class to represent numeric quantities with asymmetric uncertainties.",
          long_description=readme,
          license=license,
          install_requires=["numpy","matplotlib"])