Overview
========

HSR represents an approach to molecular similarity assessment, 
leveraging a multidimensional array to encapsulate both spatial and atomic features of molecules.
The method is grounded in a robust and deterministic process, ensuring precision and consistency in similarity comparisons.

Initial Data Representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Molecules are represented in an N-dimensional array, referred as hypershapes, where the first three dimensions correspond to 3D spatial coordinates (:func:`molecule_to_ndarray <hsr.pre_processing.molecule_to_ndarray>`).

- Additional features are integrated, enhancing the molecular description. In the default setting (:mod:`Utils <hsr.utils>`), these include:

    - The square root of the **proton** count.
    - The square root of number of **neutrons** and the number of neutrons of the most common isotope.
    - Formal charge (**electrons** information).

Principal Component Analysis (PCA) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The method applies PCA to the N-dimensional (6D in default mode) molecular representation, extracting principal components of the hyper shape (:func:`compute_pca_using_covariance <hsr.pca_transform.compute_pca_using_covariance>`).

- Orientation of eigenvectors is determined rigorously. The sign of each eigenvector is set based on the maximum projection (PCA score) of the data onto that eigenvector. This ensures a deterministic and unambiguous assignment of orientation (:func:`adjust_eigenvector_signs <hsr.pca_transform.adjust_eigenvector_signs>`).

Fingerprint Construction
~~~~~~~~~~~~~~~~~~~~~~~~
- Post-PCA, the method constructs a molecular fingerprint (:mod:`Fingerprint <hsr.fingerprint>`). This involves selecting reference points corresponding to each principal component and the geometric center of the molecule.
- The distance of each reference point from the center can be adjusted. By default, it is set to the maximum coordinate value in that dimension.
- For each reference point, distances to all atoms are calculated, resulting in a set of distance distributions.
- From each distribution, three statistical moments are computed: mean, standard deviation, and skewness. These values are collected into a list, forming the comprehensive fingerprint of the molecule.

Similarity Measurement
~~~~~~~~~~~~~~~~~~~~~~
- Molecular similarity is quantified using the normalized inverse Manhattan distance between the fingerprints of two molecules (:mod:`Similarity <hsr.similarity>`). This metric provides a straightforward yet effective measure of similarity, capturing both spatial and feature-based nuances.


Examples
~~~~~~~~

The HSR method can be used to compute the similarity between two RDKit molecules:

.. code-block:: bash

    hsr -s mol1.sdf mol2.sdf

This command will output the similarity score between the two molecules. HSR can process ``.mol``, ``.mol2``, ``.pdb``, ``.xyz``, and ``.sdf`` files.

The Manhattan distance between the fingerprints of the two molecules can also be computed:

.. code-block:: bash

    hsr -d mol1.sdf mol2.sdf

This command will output the Manhattan distance between the two molecules.

And also we can inspect the fingerprints of the molecules:

.. code-block:: bash

    hsr -f mol1.sdf mol2.sdf

This command will output the fingerprints of the two molecules. This command can take multiple molecules as input.


Three optional flags can be used to modify the behavior of the HSR tool:

- ``-chirality``: Enable chirality detection. This flag is set to ``False`` by default as chirality introduces additional complexity and potential reliability issues. For more detailed information on this aspect, please refer to our publication (TODO: add reference. Publishing in progress!).
- ``-removeH``: Remove hydrogen atoms from the molecule. This flag is set to ``False`` by default.
- ``-features FEATURES``: Possibility to choose the features to be used in the fingerprint generation. Available features are: DEFAULT_FEATURES: the default 6D representation, None: only spatial coordinates, and PROTON_FEATURES: 4D representattion of spatial coordinates and proton number. This flag is set to ``DEFAULT_FEATURES`` by default.  


Adding New Features
~~~~~~~~~~~~~~~~~~~

The HSR tool comes with its default features, but users have the flexibility to define new ones for their specific needs. 
New features must be capable of extracting or adding a property to each atom, optionally scaled as desired.

To add new features, simply define a dictionary with the new feature name as the key and a list of functions as the value.

.. code-block:: python

    EXAMPLE = {
        'new_feature_1': [extract_new_feature_1]
        'new_feature_2': [extract_new_feature_2]
        ...
    }

For comparison, here is the dictionary of the default features:

.. code-block:: python

    DEFAULT_FEATURES = {
        'protons': [extract_proton_number],
        'delta_neutrons': [extract_neutron_difference_from_common_isotope],
        'formal_charges': [extract_formal_charge]
    }

For detailed insights into the implementation and management of these features and the relative functions within HSR, refer to the :mod:`Utils <hsr.utils>` module.


Disclaimer
~~~~~~~~~~

Introducing chirality into the similarity measurement process can make the method less reliable.
For example, the PCA process can reduce the dimensionality of a molecule and hence, some of the eigenvectors's orientation will not be consistently assigned.
In this case the program will issue the following warning:

.. code-block:: python

    "WARNING: Chirality may not be consistent. {original_eigenvectors_number-len(significant_indices)} vectors have arbitrary signs."


Another case where chirality can prove unreliable is when comparing molecules with differing dimensionality, such as a different number of principal components. 
An example of this might be comparing similar 3D molecules where one has charges and the other is neutral.
In such cases, the addition of chirality detection may further reduce the similarity score. 
For detailed explanations, please refer to our publication (TODO: add reference- Publishing in progress!).

We recommend enabling chirality detection only in scenarios where molecules are unlikely to be described 
by different numbers of dimensions. However, it's important to note that this probability is hard to be 
completely eliminated, as some molecules might be planar, leading to dimensionality reduction after PCA.
Therefore, if chirality is set to `True` and the dimensionality of the two molecules being compared differs, 
the method will issue a warning as follows:

.. code-block:: python

    "WARNING: Comparison between molecules of different dimensionality: {dimensionality1} and {dimensionality2}"
                   "The similarity score may not be accurate!"


**IMPORTANT NOTE:**

   When the `chirality` parameter is set to `True`, both the :func:`compute_pca_using_covariance` and :func:`generate_fingerprint_from_molecule` functions return an additional value – the dimensionality of the molecule. This change in return values is crucial to note, especially when these methods are used in a new script.

   The :func:`compute_similarity` function is designed to handle these additional return values correctly. It will process the dimensionality information and issue a warning if there is a mismatch in dimensionality between the two molecules being compared. This is particularly important because a difference in dimensionality can significantly impact the accuracy of the similarity score.

   If you are using :func:`compute_pca_using_covariance` or :func:`generate_fingerprint_from_molecule` directly in your code, be prepared to handle an additional return value (the dimensionality) when `chirality` is `True`. This is especially relevant if you are integrating these functions into a larger workflow or using them in conjunction with other methods.

   For example, if you are performing PCA transformation step-by-step, you should modify your code to accommodate the additional dimensionality information. Similarly, when generating fingerprints, ensure that your code can handle the extra return value without errors.

   This change in the return structure is a direct consequence of enabling chirality detection, which adds a layer of complexity to the analysis but can provide more nuanced insights, especially for chiral molecules.
