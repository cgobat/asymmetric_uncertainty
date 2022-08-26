from setuptools import setup

with open("./README.md","r") as f:
    readme = f.read()

with open("./LICENSE","r") as f:
    license = f.read()

if __name__ == '__main__':
    setup(name="asymmetric_uncertainty",
          version="0.2.0",
          author="Caden Gobat",
          author_email="<cgobat@gwu.edu>",
          url="https://github.com/cgobat/asymmetric_uncertainty",
          description="A package that implements a class to represent numeric quantities with asymmetric uncertainties.",
          long_description=readme,
          license=license,
          install_requires=["numpy","matplotlib"])