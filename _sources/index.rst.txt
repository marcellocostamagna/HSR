.. HSR documentation master file, created by
   sphinx-quickstart on Mon Nov 20 14:29:53 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to HSR's documentation!
==================================

Hypershape Recognition (HSR): a Generalised Framework for Moment-Based Molecular Similarity
-------------------------------------------------------------------------------------------

HSR is a versatile, moment-based similarity measure tailored for three-dimensional (3D) chemical representations 
annotated with atomic features. It enhances the robustness and versatility of the Ultrafast Shape Recognition (USR) 
method by incorporating multidimensional features for each atom, such as protons, neutrons, and formal charges.

Getting Started
---------------

Installing HSR
~~~~~~~~~~~~~~~~

You can install HSR using either pip or conda:

.. code-block:: bash

    pip install hsr

or 

.. code-block:: bash 

    conda install hsr -c conda-forge

Build from source
--------------------

Clone the `HSR repository <https://github.com/marcellocostamagna/HSR>`_ on your machine. Move inside it and create the conda environment:


.. code-block:: bash

    conda env create -f environment.yml
    conda activate HSR_devel

Verify the correct creation of the environment by running:

.. code-block:: bash

    pytest

To use HSR from CLI run:

.. code-block:: bash

    python -m hsr.hsr_cli -h

If HSR is installed with pip or conda, the above command is replaced by the simple use of ``hsr``.

Basic Usage
------------

Run the folowing command to get help in using HSR from CLI:

.. code-block:: bash

    hsr -h

For more deatails of HSR's methodology check our :doc:`overview <detailed_overview>`.

Licensing
---------

HSR is licensed under the GNU Affero General Public License Version 3, 19 November 2007. For more details, 
see the LICENSE file in the `source code repository <https://github.com/marcellocostamagna/HSR>`_ or visit `GNU AGPL v3 License <https://www.gnu.org/licenses/agpl-3.0.html>`_.

Citing HSR
----------

If you use HSR in your research, please cite it as follows:

[TODO: Add citation]

Contributing to HSR
----------------------

We welcome contributions to HSR! If you're interested in helping, 
please read our :doc:`Contributing Guidelines <CONTRIBUTING>` for information on how to get started.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   modules
   detailed_overview
   CONTRIBUTING

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
