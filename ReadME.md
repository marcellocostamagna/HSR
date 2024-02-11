# Hyper Shape recognition (HSR): A Versatile Similarity Measure for Chemical Systems

 HSR is a versatile tool for measuring similarity across a wide range of chemical systems. It enhances the robustness and versatility of the Ultrafast Shape Recognition (USR) method by incorporating multidimensional features for each atom, such as protons, neutrons, and formal charges.

## Getting Started

### Installing ND_sim

You can install HSR using either pip or conda:

```bash
pip install hsr
```
or 

```bash
conda install hsr -c conda-forge
```

### Basic Usage

Using HSR is straightforward:

```python
import hsr

mol1 = load_molecules_from_sdf('mol1.sdf')
mol2 = load_molecules_from_sdf('mol2.sdf')

similarity = compute_similarity(mol1, mol2)

```
For a detailed overview of ND_sim's methodology check our [documentation](https://marcellocostamagna.github.io/ND_sim/).


### Licensing

HSR is licensed under the [GNU Affero General Public License Version 3](https://www.gnu.org/licenses/agpl-3.0.html). For more details, see the LICENSE file.

### Citing ND_sim

If you use HSR in your research, please cite it as follows:

[TODO: Add citation]

### Contributing

Contributions to HSR are welcome! Please read our [Contributing Guidelines](https://marcellocostamagna.github.io/ND_sim/CONTRIBUTING.html) for information on how to get started.


