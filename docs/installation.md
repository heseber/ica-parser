# Installation instructions

## Installation of the icaparser package

It is recommended to create a new virtual environment with Python >= 3.9 and to
install the icaparser package in that environment. Activate the environment and
run:

```sh
pip install "git+https://github.com/Bayer-Group/ica-parser.git"
```

If you want to install a particular development branch, use

```sh
pip install "git+https://github.com/Bayer-Group/ica-parser.git@BRANCHNAME"
```

If you use Jupyter notebooks, the virtual environment should be added as a new
Jupyter kernel. See [Using Virtual Environments in Jupyter Notebook and Python -
Parametric Thoughts](https://janakiev.com/blog/jupyter-virtual-envs/) how to do
that.

## Installation of ipywidgets

Required for progress bars in Jupyter. Please refer to the Jupyter or JupyterLab
documentation how to install the widgets. For example:

```sh
conda install jupyter # if not installed yet
conda install jupyterlab_widgets
jupyter labextension install jupyter-matplotlib
jupyter lab build
exit
```

**→ Restart Jupyter**
