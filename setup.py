from setuptools import setup

setup(  name="Kea",
        version="0.1",
        description="An object orientated histogram for use with BPASS",
        author="Max Briel",
        author_email="max.briel@auckland.ac.nz",
        packages=['kea'],
        zip_safe=False,
        install_requires=["numpy", "matplotlib"]
        )
