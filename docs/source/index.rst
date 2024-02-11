.. HSR documentation master file, created by
   sphinx-quickstart on Mon Nov 20 14:29:53 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to HSR's documentation!
==================================

Hyper Shape recognition (HSR): a generalised moment-based molecular similarity framework
-----------------------------------------------------------------------------------------

HSR is a versatile similarity measure tool applicable to a wide range of chemical systems. 
It is developed to enable fast 3D similarity measures across diverse chemicals and molecules. 
Based on the Ultrafast Shape Recognition (USR) method, HSR enhances robustness and versatility 
compared to similar previous methods.

In addition to considering the 3D coordinates of molecules, HSR incorporates additional 
features for each atom, making it a multidimensional similarity approach. While users 
have the flexibility to modify these features, the default configuration utilizes information 
based on the protons, neutrons, and formal charges of each atom.

Getting Started
---------------

Installing HSR
~~~~~~~~~~~~~~~~

HSR can be easily installed using pip or conda. Choose the method that best suits your environment:

Using pip:

.. code-block:: bash 

    pip install hsr

Using conda:

.. code-block:: bash

    conda install hsr -c conda-forge


Understanding HSR
--------------------

Although the use of HSR can be pretty straightforward:


.. code-block:: bash

    hsr -s mol1.sdf mol2.sdf

Some of the logic behind can be complex at first.
To learn more about how HSR works and its underlying methodology, see the :doc:`overview of the method <detailed_overview>`.

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
