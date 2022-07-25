from setuptools import setup

if __name__ == '__main__':
    setup(name="asymmetric_uncertainty",
          version="0.0.2",
          author="Caden Gobat",
          author_email="<cgobat@gwu.edu>",
          url="https://github.com/cgobat/asymmetric_uncertainty",
          description="A package that implements a class to represent numeric quantities with asymmetric uncertainties.",
          install_requires=["numpy","matplotlib"])