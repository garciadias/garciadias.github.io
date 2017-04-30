.. KMeans for stars documentation master file, created by
   sphinx-quickstart on Wed Sep 28 11:27:39 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
  :maxdepth: 2


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. contents::


================
KMeans for Stars
================

KMeans for stars is a python tool kit developed to apply *K*-means `(MacQueen, 
1967) <http://www-m9.ma.tum.de/foswiki/pub/WS2010/CombOptSem/kMeans.pdf>`_ 
clustering in astrophysics. This tools ware previously used in: 

- `Sánchez Almeida et al. (2016) <http://arxiv.org/abs/1601.01631>`_ - Search for Extremely Metal-poor Galaxies in the Sloan Digital Sky Survey (II): high electron temperature objects

- `Ordovás-Pascual & Sánchez Almeida (2014) <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2014A%26A...565A..53O&db_key=AST&link_type=PREPRINT>`_ - A fast version of the k-means classification algorithm for astronomical applications

- `Sánchez Almeida & Allende Prieto (2013) <http://arxiv.org/abs/1211.5321>`_ - Automated unsupervised classification of the Sloan Digital Sky Survey stellar spectra using k-means clustering

- `Morales-Luis et al. (2011) <http://arxiv.org/abs/1109.0235>`_ - Systematic search for extremely metal poor galaxies in the Sloan Digital Sky Survey

- `Sánchez Almeida et. al (2010) <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2010ApJ...714..487S&db_key=AST&link_type=PREPRINT>`_ - Automatic Unsupervised Classification of All Sloan Digital Sky Survey Data Release 7 Galaxy Spectra

- `Sánchez Almeida et. al (2009) <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2009ApJ...698.1497S&db_key=AST&link_type=PREPRINT>`_ -  Search for Blue Compact Dwarf Galaxies During Quiescence. II. Metallicities of Gas and Stars, Ages, and Star Formation Rates

- `Sánchez Almeida & Lites (2000) <http://iopscience.iop.org/article/10.1086/308603/pdf>`_ - Physical Properties of the Solar Magnetic Photosphere under the MISMA Hypothesis. II. Network and Internetwork Fields at the Disk Center

.. math::

   (a + b)^2 = a^2 + 2ab + b^2

   (\alpha - b)^2 = \alpha^2 - 2 \alpha b + b^2

Introduction
============

.. figure:: Groups_MH_Teff_cont_fill.png
   :scale: 30 %
   :alt: alternate text
   :align: center
   :figclass: align-left

   **Contour diagrams in the Teff x MH space. Different colors are used to distinguish classes, there is any touching borders with the same color. Each class is represented by five color shades, from dark to light, the shades enclose 15%, 30%, 45% and** :math:`1\sigma`. **The groups are parted in four plots minimizing superposition of classes. Panel a shows groups 0, 1, 4 and 5, panel** b **groups 7, panel** c **groups 2 and 3, and panel** d **shows group 6. Except of panel** b **, which has its classes identified as in the legend box, each class is tagged with a floating label in the form G** X **C** XX **, C referring to class and G to its group.**

Documenting objects
===================

One of Sphinx' main objectives is easy documentation of :dfn:`objects` (in a
very general sense) in any :dfn:`domain`.  A domain is a collection of object
types that belong together, complete with markup to create and reference
descriptions of these objects.

The most prominent domain is the Python domain.  To e.g. document the Python
built-in function ``enumerate()``, you would add this to one of your source
files::

   .. py:function:: enumerate(sequence[, start=0])

      Return an iterator that yields tuples of an index and an item of the
      *sequence*. (And so on.)

This is rendered like this:

.. py:function:: enumerate(sequence[, start=0])

   Return an iterator that yields tuples of an index and an item of the
   *sequence*. (And so on.)

The argument of the directive is the :dfn:`signature` of the object you
describe, the content is the documentation for it.  Multiple signatures can be
given, each in its own line.

The Python domain also happens to be the default domain, so you don't need to
prefix the markup with the domain name

Test code input
===============

Run::

	from KMeans impor assign_array
	assign_array(Centres, Your_data, i_worker, output)

.. literalinclude:: ../KMeans_garcia-dias.py
    :language: python
    :lines: 8-33

Printing some interactive coding:
---------------------------------

>>> 1+1
2

Some text that requires a footnote [#f1]_ .

.. rubric:: Footnotes

.. [#f1] Text of the first footnote.
