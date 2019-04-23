import setuptools

setuptools.setup(
    name='ambrose',
    version='2.0',
    author='mszegedy',
    author_email='mszegedy2@gmail.com',
    description='AMBER interoperability for PyRosetta.',
    url='https://github.com/RosettaCommons/tools/tree/master/AmbRose',
    packages=setuptools.find_packages(exclude=['legacy']),
    classifiers=[
        'Programming Language :: Python :: 3'
    ]
)
