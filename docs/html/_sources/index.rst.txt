======================
Drinfeld Modular Forms
======================

This is a SageMath package for computing with Drinfeld modular forms and their
expansion at infinity.

This work is licensed under a `GNU General Public License <https://www.gnu.org/licenses/>`_.

Installation
============

This package has been tested on SageMath version 9.8 and higher. It is
not guaranteed to work on previous versions.

* Install from PyPI

The easiest way to install this package is via the Python package index.
You simply have to run SageMath first and then type the following
command::

   sage: pip install drinfeld-modular-forms

* Install from source code

You can also install this package by cloning the source code from the
`Github repo <https://github.com/DavidAyotte/drinfeld_modular_forms>`_.

Next, you have to run ``make install`` inside the project's folder. You
can also run the following command::

   sage -pip install --upgrade --no-index -v .

If any updates are made to the package, then you have to pull them from
the repo and run the above command again.

Usage
=====

After running SageMath, you can import the functionalities of this
package by typing the following command::

   sage: from drinfeld_modular_forms import *

Table of Contents
=================

.. toctree::
   :maxdepth: 2

   background
   ring
   element
   expansions
   goss_polynomials
   drinfeld_modules


Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
