# Building the LIGGGHTS documentation with Sphinx

This document describes how to build the [Sphinx-based](http://www.sphinx-doc.org) documentation of LIGGGHTS.

## Prerequisites

### Python 3

The [LAMMPS Documentation Utilities](https://github.com/ParticulateFlow/lammps-doc-utils), which are discussed below, require Python 3. Thus, the LIGGGHTS documentation requires [Python 3](https://docs.python.org/3/). Apart from Python 3 itself, [pip](https://pypi.python.org/pypi/pip) is required to install some of the prerequisites.

On an Ubuntu system run the following command to install *pip* and all its dependencies.

    sudo apt-get install python3-pip


### Sphinx

[Sphinx](http://www.sphinx-doc.org) is obviously required. Furthermore, the package [sphinxcontrib.images](https://github.com/spinus/sphinxcontrib-images) is used by throughout the documentation to deal with images. Both, Sphinx itself and the additional package can be installed using *pip*

    pip3 install sphinx
    pip3 install sphinxcontrib-images

Note that `pip` defaults on most systems to *pip* of Python 2, thus `pip3` is explicitly called.


### LAMMPS Documentation Utilities

The [Sphinx-based](http://www.sphinx-doc.org) documentation of LIGGGHTS is based on the [LAMMPS Documentation Utilities](https://github.com/ParticulateFlow/lammps-doc-utils).
These contain the tool *txt2rst*, which converts the files in the plain-text, LAMMPS documentation format to Sphinx ReStructured Text.

First clone the [repository](https://github.com/ParticulateFlow/lammps-doc-utils), and then run the installer

    git clone https://github.com/ParticulateFlow/lammps-doc-utils
    cd lammps-doc-utils
    python3 setup.py install

Note that `python` defaults to Python 2 on many systems, that's why `python3` is explicitly called in the last line.


## Building the documentation

Building the actual documentation is fairly easy. Just change into the `doc` folder of your LIGGGHTS installation and run the command `make` with the appropriate target (e.g. html or latexpdf). Running `make` without a target will print a list of all available targets.


    cd $CFDEM_LIGGGHTS_INST_DIR
    cd doc
    make html

The default settings in the `Makefile` will cause Sphinx to create the documentation at
`$CFDEM_LIGGGHTS_INST_DIR/doc/_build/<target>`.

Thus, the path to the HTML documentation will be
`$CFDEM_LIGGGHTS_INST_DIR/doc/_build/html/Manual.html`


### Trouble-shooting

#### Choosing the correct Sphinx version

If there are multiple versions of Sphinx present on your system, e.g. Sphinx for Python 2 and Python 3, the command `sphinx-build` might default to the wrong version. In this case you can explicitly state the full path to the command `sphinx-build` of Python 3 at the top of the `Makefile`, e.g.:

    SPHINXBUILD   = /usr/share/sphinx/scripts/python3/sphinx-build
