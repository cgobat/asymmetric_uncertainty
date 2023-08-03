from setuptools import setup, find_packages

with open("./README.md","r") as f:
    readme = f.read()

if __name__ == '__main__':
    setup(name="asymmetric_uncertainty",
          version="0.2.2",
          author="Caden Gobat",
          author_email="cgobat@gwu.edu",
          packages=find_packages(),
          url="https://github.com/cgobat/asymmetric_uncertainty",
          description="A package that implements a class to represent numeric quantities with asymmetric uncertainties.",
          long_description=readme,
          download_url="https://github.com/cgobat/asymmetric_uncertainty/archive/main.zip",
          license="GPL-3.0",
          install_requires=["numpy", "matplotlib"])